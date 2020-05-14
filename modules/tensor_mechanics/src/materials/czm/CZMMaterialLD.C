//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CZMMaterialLD.h"
#include "RankFourTensor.h"
#include "RotationMatrix.h"

registerMooseObject("MooseApp", CZMMaterialLD);

template <>
InputParameters
validParams<CZMMaterialLD>()
{
  InputParameters params = InterfaceMaterial::validParams();
  params.addClassDescription("Base class for cohesive zone material models for rate depndent");
  params.addRequiredCoupledVar("displacements",
                               "The string of displacements suitable for the problem statement");
  params.addClassDescription("Base czm material for large deformation");
  params.addParam<bool>("large_kinematics", true, "Use large displacement kinematics.");
  params.addParam<Real>("E", 1e3, "opening stiffness");
  params.addParam<Real>("G", 1e2, "shear modulus");
  params.addParam<bool>(
      "rotate_traction_back", false, "Formulate residual in reference configuration");
  params.addParam<bool>("use_area_change", false, "account for interface area changes");
  return params;
}

CZMMaterialLD::CZMMaterialLD(const InputParameters & parameters)
  : InterfaceMaterial(parameters),
    _normals(_assembly.normals()),
    _ndisp(coupledComponents("displacements")),
    _disp(3),
    _disp_neighbor(3),
    _disp_vars(3),
    _disp_old(3),
    _disp_neighbor_old(3),
    _displacement_jump_global(declareProperty<RealVectorValue>("displacement_jump_global")),
    _displacement_jump_global_old(declareProperty<RealVectorValue>("displacement_jump_global_old")),
    _displacement_jump(declareProperty<RealVectorValue>("displacement_jump")),
    _displacement_jump_old(getMaterialPropertyOld<RealVectorValue>("displacement_jump")),
    _displacement_jump_global_inc(declareProperty<RealVectorValue>("displacement_jump_global_inc")),
    _displacement_jump_inc(declareProperty<RealVectorValue>("displacement_jump_inc")),

    _traction_global(declareProperty<RealVectorValue>("traction_global")),
    _traction_global_inc(declareProperty<RealVectorValue>("traction_global_inc")),
    _traction(declareProperty<RealVectorValue>("traction")),
    _traction_inc(declareProperty<RealVectorValue>("traction_inc")),
    _traction_global_old(getMaterialPropertyOld<RealVectorValue>("traction_global")),
    _traction_old(getMaterialPropertyOld<RealVectorValue>("traction")),
    _traction_derivatives_global(
        declareProperty<std::vector<RankTwoTensor>>("traction_derivatives_global")),
    _traction_derivatives_global_neighbor(
        declareProperty<std::vector<RankTwoTensor>>("traction_derivatives_global_neighbor")),
    _traction_derivatives(declareProperty<RankTwoTensor>("traction_derivatives")),

    _F_avg(declareProperty<RankTwoTensor>("F_avg")),
    _DF_avg(declareProperty<RankTwoTensor>("DF_avg")),
    _F_avg_old(getMaterialPropertyOld<RankTwoTensor>("F_avg")),
    _R_avg(declareProperty<RankTwoTensor>("R_avg")),
    _R_avg_old(getMaterialPropertyOld<RankTwoTensor>("R_avg")),
    _U_avg(declareProperty<RankTwoTensor>("U_avg")),
    _DR_avg(declareProperty<RankTwoTensor>("DR_avg")),
    _DL_avg(declareProperty<RankTwoTensor>("DL_avg")),
    _DL_avg_old(getMaterialPropertyOld<RankTwoTensor>("DL_avg")),
    _n_avg(declareProperty<RealVectorValue>("n_avg")),
    _dsdotdS_avg(declareProperty<Real>("dsdotdS_avg")),
    _dadA_avg(declareProperty<Real>("dadA_avg")),

    _ld(getParam<bool>("large_kinematics")),
    _use_area_change(getParam<bool>("use_area_change"))
{ // Enforce consistency
  if (_ndisp != _mesh.dimension())
    paramError("displacements", "Number of displacements must match problem dimension.");
}

void
CZMMaterialLD::initialSetup()
{
  _K(0, 0) = getParam<Real>("E");
  _K(1, 1) = getParam<Real>("G");
  _K(2, 2) = getParam<Real>("G");
  // initializing the displacement vectors
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    _disp[i] = &coupledValue("displacements", i);
    _disp_neighbor[i] = &coupledNeighborValue("displacements", i);
    _disp_vars[i] = getVar("displacements", i);
    _disp_old[i] = &coupledValueOld("displacements", i);
    _disp_neighbor_old[i] = &coupledNeighborValueOld("displacements", i);

    _grad_disp.push_back(&coupledGradient("displacements", i));
    _grad_disp_neighbor.push_back(&coupledNeighborGradient("displacements", i));
    _phi.push_back(&(_disp_vars[i]->phiFace()));
    _phi_neighbor.push_back(&(_disp_vars[i]->phiFaceNeighbor()));
    _grad_phi.push_back(&(_disp_vars[i]->gradPhiFace()));
    _grad_phi_neighbor.push_back(&(_disp_vars[i]->gradPhiFaceNeighbor()));
  }

  // All others zero (so this will work naturally for 2D and 1D problems)
  for (unsigned int i = _ndisp; i < 3; i++)
  {
    _disp[i] = &_zero;
    _disp_neighbor[i] = &_zero;
    _disp_old[i] = &_zero;
    _disp_neighbor_old[i] = &_zero;
    _grad_disp.push_back(&_grad_zero);
    _grad_disp_neighbor.push_back(&_grad_zero);
    _phi.push_back(&_phi_zero);
    _phi_neighbor.push_back(&_phi_zero);
    _grad_phi.push_back(&_grad_phi_zero);
    _grad_phi_neighbor.push_back(&_grad_phi_zero);
  }
}

void
CZMMaterialLD::initQpStatefulProperties()
{

  _traction_global[_qp] = 0;
  _traction[_qp] = 0;
  _displacement_jump[_qp] = 0;
  _F_avg[_qp] = RankTwoTensor::Identity();
  _R_avg[_qp] = RankTwoTensor::Identity();
  _DL_avg[_qp] = RankTwoTensor::initNone;
}

void
CZMMaterialLD::computeQpProperties()
{
  const unsigned int n_sf = _disp_vars[0]->gradPhi().size();
  _djump_global_duk.resize(n_sf);
  _djump_global_duk_neighbor.resize(n_sf);
  _dFmn_duk.resize(n_sf);
  _dFmn_neighbor_duk.resize(n_sf);

  for (unsigned int sf = 0; sf < n_sf; sf++)
  {
    for (unsigned int i = 0; i < 3; i++)
    {
      _djump_global_duk[sf](i, i) = _disp_vars[i]->phiFace()[sf][_qp];
      _djump_global_duk_neighbor[sf](i, i) = _disp_vars[i]->phiFaceNeighbor()[sf][_qp];
    }
    if (_qp == 0)
    {
      std::cout << " sf" << sf << std::endl;
      std::cout << " phi \n" << _djump_global_duk[sf];
      std::cout << " phineighbor \n" << _djump_global_duk_neighbor[sf];
      std::cout << "\n";
    }
  }
  std::cout << "\n\n\n\n";
  for (unsigned int r = 0; r < n_sf; r++)
    for (unsigned int i = 0; i < 3; i++)
    {
      for (unsigned int j = 0; j < 3; j++)
      {
        _dFmn_duk[r](i, j, i) = 0.5 * (*_grad_phi[i])[r][_qp](j);
        _dFmn_neighbor_duk[r](i, j, i) = (*_grad_phi_neighbor[i])[r][_qp](j);
      }
    }

  // update the displacement jump
  for (unsigned int i = 0; i < 3; i++)
  {
    _displacement_jump_global[_qp](i) = (*_disp_neighbor[i])[_qp] - (*_disp[i])[_qp];
    _displacement_jump_global_old[_qp](i) = (*_disp_neighbor_old[i])[_qp] - (*_disp_old[i])[_qp];
  }
  _displacement_jump_global_inc[_qp] =
      (_displacement_jump_global[_qp] - _displacement_jump_global_old[_qp]);

  const RealTensorValue Q0 = RotationMatrix::rotVec1ToVec2(_normals[_qp], RealVectorValue(1, 0, 0));
  // RankTwoTensor DR, DF, DRF, DUF, Rdot;

  _n_avg[_qp] = _normals[_qp];
  _dsdotdS_avg[_qp] = 0.;
  _dadA_avg[_qp] = 1.;
  _DR_avg[_qp] = RankTwoTensor::initNone;
  _R_avg[_qp] = RankTwoTensor::Identity();

  if (_ld)
  {
    computeFandL();
    PolarDecomposition(_F_avg[_qp], _R_avg[_qp], _U_avg[_qp]);
    _DF_avg[_qp] = _F_avg[_qp] * _F_avg_old[_qp].inverse();
    _DR_avg[_qp] = _R_avg[_qp] - _R_avg_old[_qp];

    if (_use_area_change)
    {
      _n_avg[_qp] = _R_avg[_qp] * _normals[_qp];
      _dsdotdS_avg[_qp] = _DL_avg[_qp].trace() - _n_avg[_qp] * (_DL_avg[_qp] * _n_avg[_qp]);
      _dadA_avg[_qp] =
          _F_avg[_qp].det() * ((_F_avg[_qp].inverse().transpose() * _normals[_qp]).norm());
    }
  }

  _displacement_jump[_qp] = Q0 * _R_avg[_qp].transpose() * _displacement_jump_global[_qp];
  _displacement_jump_inc[_qp] = Q0 * (_DR_avg[_qp].transpose() * _displacement_jump_global[_qp] +
                                      _R_avg[_qp].transpose() * _displacement_jump_global_inc[_qp]);

  computeTractionAndDerivatives();

  // assemble PK1 traction
  _traction_global_inc[_qp] = _DR_avg[_qp] * Q0.transpose() * _traction[_qp] +
                              _R_avg[_qp] * Q0.transpose() * _traction_inc[_qp];

  if (_ld && _use_area_change)
  {
    _traction_global_inc[_qp] += _dsdotdS_avg[_qp] * _R_avg[_qp] * Q0.transpose() * _traction[_qp];
    _traction_global_inc[_qp] *= _dadA_avg[_qp];
  }

  _traction_global[_qp] = _traction_global_old[_qp] + _traction_global_inc[_qp];

  // let's start with derivatives
  RankFourTensor dRdF;
  Real a = 1;
  RankTwoTensor B = RankTwoTensor::initNone;
  RankTwoTensor C = RankTwoTensor::initNone;
  RankTwoTensor D = _R_avg[_qp] * Q0.transpose();

  std::vector<RealVectorValue> da_du(n_sf);
  std::vector<RankThreeTensor> dB_du(n_sf);
  std::vector<RankThreeTensor> dC_du(n_sf);
  std::vector<RankThreeTensor> dD_du(n_sf);

  std::vector<RealVectorValue> da_du_neighbor(n_sf);
  std::vector<RankThreeTensor> dB_du_neighbor(n_sf);
  std::vector<RankThreeTensor> dC_du_neighbor(n_sf);
  std::vector<RankThreeTensor> dD_du_neighbor(n_sf);

  dRdF.zero();
  if (_ld)
  {
    dRdF = computedRdF(_R_avg[_qp], _U_avg[_qp]);
    C = _DR_avg[_qp] * Q0.transpose();
    if (_use_area_change)
    {
      a = _dadA_avg[_qp];
      B = _dsdotdS_avg[_qp];
    }
  }

  _traction_derivatives_global[_qp].resize(n_sf);
  _traction_derivatives_global_neighbor[_qp].resize(n_sf);

  // std::cout << "a: " << a << std::endl;
  // std::cout << "B: " << B << std::endl;
  // std::cout << "C: " << C << std::endl;
  // std::cout << "D: " << D << std::endl;

  // std::cout << "K " << _R_avg[_qp] * Q0.transpose() * _K * Q0 * _R_avg[_qp].transpose();

  const RankTwoTensor RpDR = (_R_avg[_qp] + _DR_avg[_qp]).transpose();
  //
  // std::cout << "computed K" << a * ((B + C + D) * _traction_derivatives[_qp] * Q0 * RpDR);

  RankTwoTensor temp;
  for (unsigned int sf = 0; sf < n_sf; sf++)
  {
    RankTwoTensor djump_du = Q0 * RpDR * _djump_global_duk[sf];
    RankTwoTensor djump_du_neighbor = Q0 * RpDR * _djump_global_duk_neighbor[sf];
    // _traction_derivatives_global_neighbor[_qp][sf] =
    //     da_du_neighbor[sf] * ((B + C) * _traction_global[_qp] + D * _traction_global_inc[_qp]);
    //
    // _traction_derivatives_global[_qp][sf] =
    //     da_du[sf] * ((B + C) * _traction_global[_qp] + D * _traction_global_inc[_qp]);
    //
    // temp = a * (RijkVk(dB_du_neighbor[sf] + dC_du_neighbor[sf], _traction_global[_qp]) +
    //             RijkVk(dC_du_neighbor[sf], _traction_global_inc[_qp]));
    // _traction_derivatives_global_neighbor[_qp][sf] += temp;
    //
    // temp = a * (RijkVk(dB_du[sf] + dC_du[sf], _traction_global[_qp]) +
    //             RijkVk(dC_du[sf], _traction_global_inc[_qp]));
    // _traction_derivatives_global[_qp][sf] += temp;

    temp = a * ((B + C + D) * (_traction_derivatives[_qp] * djump_du_neighbor));
    _traction_derivatives_global_neighbor[_qp][sf] += temp;

    temp = a * ((B + C + D) * (_traction_derivatives[_qp] * djump_du));
    _traction_derivatives_global[_qp][sf] += temp;
  }

  // RankFourTensor dRdF = RankFourTensor::initNone;
  // RankTwoTensor djumpelem_dFmn = Q0;
  // RankTwoTensor djumpelem_phi = Q0;
  // if (_ld) {
  //   dRdF = computedRdF(_R_avg[_qp], _U_avg[_qp]);
  //   djumpelem_dFmn *= (_DR_avg[_qp] + _R_avg[_qp]).transpose();
  //
  //   // compute djumpelem_phi term
  //   djumpelem_phi *= (-dRdF * (_displacement_jump_global[_qp] +
  //                              _displacement_jump_global_inc[_qp]))
  //       _traction_derivatives[_qp] =
  //           _R_avg[_qp] * Q0.transpose() * _traction_derivatives[_qp]
  // *
  //           djumpelem_djumpglobal;
  // }
  // }
  // }
}

void
CZMMaterialLD::computeFandL()
{

  RankTwoTensor F =
      (RankTwoTensor::Identity() +
       RankTwoTensor((*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp], (*_grad_disp[2])[_qp]));
  RankTwoTensor F_neighbor =
      (RankTwoTensor::Identity() + RankTwoTensor((*_grad_disp_neighbor[0])[_qp],
                                                 (*_grad_disp_neighbor[1])[_qp],
                                                 (*_grad_disp_neighbor[2])[_qp]));

  _F_avg[_qp] = 0.5 * (F + F_neighbor);
  _DL_avg[_qp] = RankTwoTensor::Identity() - _F_avg_old[_qp] * _F_avg[_qp].inverse();
}

void
CZMMaterialLD::PolarDecomposition(const RankTwoTensor & F,
                                  RankTwoTensor & R,
                                  RankTwoTensor & U) const
{

  // copy paste from tnesor mechanics
  std::vector<Real> e_value(3);
  RankTwoTensor e_vector, N1, N2, N3;

  RankTwoTensor C = F.transpose() * F;
  C.symmetricEigenvaluesEigenvectors(e_value, e_vector);

  const Real lambda1 = std::sqrt(e_value[0]);
  const Real lambda2 = std::sqrt(e_value[1]);
  const Real lambda3 = std::sqrt(e_value[2]);

  N1.vectorOuterProduct(e_vector.column(0), e_vector.column(0));
  N2.vectorOuterProduct(e_vector.column(1), e_vector.column(1));
  N3.vectorOuterProduct(e_vector.column(2), e_vector.column(2));

  U = N1 * lambda1 + N2 * lambda2 + N3 * lambda3;
  R = F * U.inverse();
}

// void CZMMaterialLD::computeRdot(const RankTwoTensor &F, const
// RankTwoTensor &L,
//                                 const RankTwoTensor &R, const
//                                 RankTwoTensor &U, RankTwoTensor &Rdot) {
//
//   RankTwoTensor Fdot = L * F;
//   RankTwoTensor Uhat = U.trace() * RankTwoTensor::Identity() - U;
//   Rdot = 1. / Uhat.det() * (R * Uhat) *
//          (R.transpose() * Fdot - Fdot.transpose() * R) *
//          Uhat.transpose();
// }

void
CZMMaterialLD::computeSdot(const RankTwoTensor & L, const RealVectorValue & n, Real & sdot)
{
  Real t1 = n * (L * n);
  sdot = L.trace() - t1;
}

void
CZMMaterialLD::computeTractionAndDerivatives()
{
  _traction_inc[_qp] = _K * _displacement_jump_inc[_qp];
  _traction[_qp] = _traction_old[_qp] + _traction_inc[_qp];

  _traction_derivatives[_qp] = _K;
}

// void CZMMaterialLD::computeTractionDerivative() {
//   _traction_derivatives[_qp] = _K;
// }

RankFourTensor
CZMMaterialLD::computedRdF(const RankTwoTensor & R, const RankTwoTensor & U)
{
  const RankTwoTensor Uhat = U.trace() * RankTwoTensor::Identity() - U;
  unsigned int k, l, m, n, p, q;
  const Real Uhat_det = Uhat.det();
  RankFourTensor dRdF;
  for (k = 0; k < 3; k++)
    for (l = 0; l < 3; l++)
      for (m = 0; m < 3; m++)
        for (n = 0; n < 3; n++)
        {
          for (p = 0; p < 3; p++)
          {
            for (q = 0; q < 3; q++)
              dRdF(k, l, m, n) +=
                  R(k, p) * (Uhat(p, q) * R(m, q) * Uhat(n, l) - Uhat(p, n) * R(m, q) * Uhat(q, l));
          }
          dRdF /= Uhat_det;
        }

  return dRdF;
}
