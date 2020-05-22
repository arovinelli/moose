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
    _dtractionglobal_djumpglobal(declareProperty<RankTwoTensor>("dtractionglobal_djumpglobal")),
    _dtractionglobal_dF(declareProperty<RankThreeTensor>("dtractionglobal_dF")),
    _dtraction_djump(declareProperty<RankTwoTensor>("dtraction_djump")),

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

  _perturbedF = RankTwoTensor::initNone;
  _perturbedU.zero();
  const Real delta = 1e-6;
  const Real tol = 1e-4;
  evaluateMP();
  const RealVectorValue actual_traction = _traction_global[_qp];
  const RankThreeTensor dt_df_an = _dtractionglobal_dF[_qp];
  const RankTwoTensor dt_du_an = _dtractionglobal_djumpglobal[_qp];
  const RankTwoTensor actual_R = _R_avg[_qp];
  const RankFourTensor dR_df_an = _dR_dF;

  RankThreeTensor dt_df_fd;
  _perturbedF = RankTwoTensor::initNone;
  _perturbedU.zero();
  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
    {
      _perturbedF = RankTwoTensor::initNone;
      _perturbedF(i, j) = delta;
      evaluateMP();
      for (unsigned int k = 0; k < 3; k++)
        dt_df_fd(k, i, j) = (_traction_global[_qp](k) - actual_traction(k)) / _perturbedF(i, j);
    }

  RankFourTensor dR_df_fd;
  _perturbedF = RankTwoTensor::initNone;
  _perturbedU.zero();
  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
    {
      _perturbedF = RankTwoTensor::initNone;
      _perturbedF(i, j) = delta;
      evaluateMP();
      for (unsigned int k = 0; k < 3; k++)
        for (unsigned int l = 0; l < 3; l++)
          dR_df_fd(k, l, i, j) = (_R_avg[_qp](k, l) - actual_R(k, l)) / _perturbedF(i, j);
    }

  _perturbedF = RankTwoTensor::initNone;
  _perturbedU.zero();
  RankTwoTensor dt_du_fd;
  for (unsigned int i = 0; i < 3; i++)
  {
    _perturbedU.zero();
    _perturbedU(i) = delta;
    evaluateMP();
    for (unsigned int j = 0; j < 3; j++)
      dt_du_fd(j, i) = (_traction_global[_qp](j) - actual_traction(j)) / _perturbedU(i);
  }

  if (_qp == 0)
  {
    RankThreeTensor derivative_difference_R3;
    RankFourTensor derivative_difference_R4;
    RankTwoTensor derivative_difference_R2 = RankTwoTensor::initNone;
    derivative_difference_R3.zero();
    bool print = false;
    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int j = 0; j < 3; j++)
        for (unsigned int k = 0; k < 3; k++)
        {
          derivative_difference_R3(i, j, k) = std::abs(dt_df_an(i, j, k) - dt_df_fd(i, j, k));
          print = print || derivative_difference_R3(i, j, k) > tol;
        }

    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int j = 0; j < 3; j++)
        for (unsigned int k = 0; k < 3; k++)
          for (unsigned int l = 0; l < 3; l++)
          {
            derivative_difference_R4(i, j, k, l) =
                std::abs(dR_df_an(i, j, k, l) - dR_df_fd(i, j, k, l));
            print = print || derivative_difference_R4(i, j, k, l) > tol;
          }

    if (print && _ld)
    {
      std::cout << "****** DERIVATIVE ERROR*******\n";
      std::cout << " _traction_global_DF \n";
      for (unsigned int i = 0; i < 3; i++)
        for (unsigned int j = 0; j < 3; j++)
          for (unsigned int k = 0; k < 3; k++)
          {
            if (derivative_difference_R3(i, j, k) > tol)
              std::cout << "component(" << i << "," << j << "," << k
                        << ")= " << derivative_difference_R3(i, j, k) << std::endl;
          }
      std::cout << "\n\n";

      std::cout << " _DR_DF \n";
      for (unsigned int i = 0; i < 3; i++)
        for (unsigned int j = 0; j < 3; j++)
          for (unsigned int k = 0; k < 3; k++)
            for (unsigned int l = 0; l < 3; l++)
            {
              if (derivative_difference_R4(i, j, k, l) > tol)
                std::cout << "component(" << i << "," << j << "," << k << "," << l
                          << ")= " << derivative_difference_R4(i, j, k, l)
                          << " AN= " << dR_df_an(i, j, k, l) << " FD= " << dR_df_fd(i, j, k, l)
                          << std::endl;
            }
      std::cout << "\n\n";
    }

    print = false;
    derivative_difference_R2.zero();
    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int j = 0; j < 3; j++)
      {
        derivative_difference_R2(i, j) = std::abs(dt_du_an(i, j) - dt_du_fd(i, j));
        print = print || derivative_difference_R2(i, j) > delta;
      }

    if (print && !_ld)
    {
      std::cout << "****** DERIVATIVE ERROR*******\n";
      std::cout << " _traction_global_du \n";
      for (unsigned int i = 0; i < 3; i++)
        for (unsigned int j = 0; j < 3; j++)
        {
          // if (derivative_difference_R2(i, j) > delta)
          std::cout << "component(" << i << "," << j << ")= " << derivative_difference_R2(i, j)
                    << " AN= " << dt_du_an(i, j) << " FD= " << dt_du_fd(i, j) << std::endl;
        }
      std::cout << "\n\n";
    }

    _perturbedF = RankTwoTensor::initNone;
    evaluateMP();
  }
}

void
CZMMaterialLD::evaluateMP()
{

  // update the displacement jump
  for (unsigned int i = 0; i < 3; i++)
  {
    _displacement_jump_global[_qp](i) = (*_disp_neighbor[i])[_qp] - (*_disp[i])[_qp];
    _displacement_jump_global_old[_qp](i) = (*_disp_neighbor_old[i])[_qp] - (*_disp_old[i])[_qp];
  }
  _displacement_jump_global[_qp] += _perturbedU;
  _displacement_jump_global_inc[_qp] =
      (_displacement_jump_global[_qp] - _displacement_jump_global_old[_qp]);

  const RealTensorValue Q0 = RotationMatrix::rotVec1ToVec2(_normals[_qp], RealVectorValue(1, 0, 0));

  _n_avg[_qp] = _normals[_qp];
  _dsdotdS_avg[_qp] = 0.;
  _dadA_avg[_qp] = 1.;
  _DR_avg[_qp] = RankTwoTensor::initNone;
  _R_avg[_qp] = RankTwoTensor::Identity();

  /*************************
  COMPUTE GLOBAL TRACTION
  *************************/

  // initialize large deformation useful terms
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

  const Real a = _dadA_avg[_qp];
  const RankTwoTensor C = _DR_avg[_qp] * Q0.transpose();
  const RankTwoTensor D = _R_avg[_qp] * Q0.transpose();
  const RankTwoTensor B = _dsdotdS_avg[_qp] * D;

  // compute local displacement jump
  _displacement_jump[_qp] = D.transpose() * _displacement_jump_global[_qp];
  _displacement_jump_inc[_qp] = C.transpose() * _displacement_jump_global[_qp] +
                                D.transpose() * _displacement_jump_global_inc[_qp];

  computeTractionAndDerivatives();

  // assemble PK1 traction
  _traction_global_inc[_qp] = C * _traction[_qp] + D * _traction_inc[_qp];

  if (_ld && _use_area_change)
  {
    _traction_global_inc[_qp] += B * _traction[_qp];
    _traction_global_inc[_qp] *= a;
  }

  _traction_global[_qp] = _traction_global_old[_qp] + _traction_global_inc[_qp];

  /*****************************************************
  COMPUTE DERIVATIVE W.R.T. THE GLOBAL DISAPLCMENT JUMP
  ******************************************************/

  // here we compute the derivative of the global traction wrt the actual disaplcment jump
  const RankTwoTensor djump_djumpglobal = (C + D).transpose();
  const RankTwoTensor dtraction_djumpglobal = _dtraction_djump[_qp] * djump_djumpglobal;
  _dtractionglobal_djumpglobal[_qp] = a * ((B + C + D) * dtraction_djumpglobal);

  /*****************************************************
  COMPUTE DERIVATIVE W.R.T. THE DEFORMATION GRADIENT
  ******************************************************/
  if (!_ld)
    _dtractionglobal_dF[_qp].zero();
  else
  {
    _dR_dF = computedRdF(_R_avg[_qp], _U_avg[_qp]);
    // const RankFourTensor dD_dF = dC_dF;
    RankFourTensor dC_dF, dC_dF_T;
    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int j = 0; j < 3; j++)
        for (unsigned int l = 0; l < 3; l++)
          for (unsigned int m = 0; m < 3; m++)
          {
            dC_dF(i, j, l, m) = 0;
            dC_dF_T(i, j, l, m) = 0;
            for (unsigned int k = 0; k < 3; k++)
            {
              dC_dF(i, j, l, m) += _dR_dF(i, k, l, m) * Q0(j, k);
              dC_dF_T(i, j, l, m) += Q0(i, k) * _dR_dF(j, k, l, m);
            }
          }
    // const RankFourTensor dD_dF = dC_dF;

    RankTwoTensor da_dF;
    RankFourTensor dB_dF;
    RealVectorValue tempV;
    RankTwoTensor R2temp;
    RankThreeTensor R3temp;
    RankFourTensor R4temp;

    if (_use_area_change)
    {
      // dadF
      const RankTwoTensor F_inv = _F_avg[_qp].inverse();
      const RankTwoTensor F_itr = F_inv.transpose();
      const RealVectorValue Fitr_N = F_itr * _normals[_qp];
      const Real Fitr_N_norm = Fitr_N.norm();
      const RankFourTensor dFinv_dF = dR2inverse(F_inv);

      for (unsigned int l = 0; l < 3; l++)
        for (unsigned int m = 0; m < 3; m++)
        {
          R2temp(l, m) = 0;
          for (unsigned int i = 0; i < 3; i++)
            for (unsigned int j = 0; j < 3; j++)
              R2temp(l, m) += Fitr_N(i) * dFinv_dF(j, i, l, m) * _normals[_qp](j);
        }
      da_dF = _F_avg[_qp].ddet() * Fitr_N_norm + _F_avg[_qp].det() / Fitr_N_norm * R2temp;

      // dBdF

      // now cpmpute the derivatitve of dtrace(L)dF
      // start with dDL_dF
      RankFourTensor dDL_dF;
      for (unsigned int i = 0; i < 3; i++)
        for (unsigned int j = 0; j < 3; j++)
          for (unsigned int l = 0; l < 3; l++)
            for (unsigned int m = 0; m < 3; m++)
            {
              dDL_dF(i, j, l, m) = 0;
              for (unsigned int k = 0; k < 3; k++)
                dDL_dF(i, j, l, m) -= _F_avg_old[_qp](i, k) * dFinv_dF(k, j, l, m);
            }

      RankTwoTensor dtraceDL_dF;
      R2temp = RankTwoTensor::Identity();
      for (unsigned int l = 0; l < 3; l++)
        for (unsigned int m = 0; m < 3; m++)
        {
          dtraceDL_dF(l, m) = 0;
          for (unsigned int i = 0; i < 3; i++)
            for (unsigned int j = 0; j < 3; j++)
              dtraceDL_dF(l, m) += R2temp(i, j) * dDL_dF(i, j, l, m);
        }

      // now compute the derivatitve of the the second part of dsdotdS_avg (n_i*DL_ij*n_j)
      for (unsigned int p = 0; p < 3; p++)
        for (unsigned int q = 0; q < 3; q++)
        {
          R2temp(p, q) = 0;
          for (unsigned int r = 0; r < 3; r++)
            for (unsigned int l = 0; l < 3; l++)
              for (unsigned int k = 0; k < 3; k++)
                R2temp(p, q) +=
                    _dR_dF(k, r, p, q) * _normals[_qp](r) * _DL_avg[_qp](k, l) * _n_avg[_qp](l) +
                    _dR_dF(l, r, p, q) * _normals[_qp](r) * _DL_avg[_qp](k, l) * _n_avg[_qp](k) +
                    dDL_dF(k, l, p, q) * _n_avg[_qp](k) * _n_avg[_qp](l);
        }

      // assemble ddsdotdS_dF*D
      for (unsigned int i = 0; i < 3; i++)
        for (unsigned int j = 0; j < 3; j++)
          for (unsigned int r = 0; r < 3; r++)
            for (unsigned int s = 0; s < 3; s++)
              R4temp(i, j, r, s) = (dtraceDL_dF(r, s) - R2temp(r, s)) * D(i, j);
      dB_dF = R4temp + _dsdotdS_avg[_qp] * dC_dF;
    }
    RankThreeTensor T1, T2, T3;

    // assmble T1
    if (_use_area_change)
    {
      tempV = (B + C) * _traction[_qp] + D * _traction_inc[_qp];
      for (unsigned int i = 0; i < 3; i++)
        for (unsigned int j = 0; j < 3; j++)
          for (unsigned int k = 0; k < 3; k++)
            T1(i, j, k) = da_dF(j, k) * tempV(i);
    }
    else
      T1.zero();

    // assemeble T2;
    tempV = _traction[_qp] + _traction_inc[_qp];
    T2 = a * RijklVj(dC_dF, tempV);
    if (_use_area_change)
      T2 += a * RijklVj(dB_dF, _traction[_qp]);

    // assemeble T3;
    R2temp = (B + C + D) * _dtraction_djump[_qp];
    tempV = _displacement_jump_global[_qp] + _displacement_jump_global_inc[_qp];
    R3temp = RijklVj(dC_dF_T, tempV);
    T3 = a * RijRjkl(R2temp, R3temp);

    _dtractionglobal_dF[_qp] = T1 + T2 + T3;
  }
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
  _F_avg[_qp] += _perturbedF;
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
  _dtraction_djump[_qp] = _K;
}

RankFourTensor
CZMMaterialLD::computedRdF(const RankTwoTensor & R, const RankTwoTensor & U)
{
  const RankTwoTensor Uhat = U.trace() * RankTwoTensor::Identity() - U;
  unsigned int k, l, m, n, p, q;
  const Real Uhat_det = Uhat.det();

  for (k = 0; k < 3; k++)
    for (l = 0; l < 3; l++)
      for (m = 0; m < 3; m++)
        for (n = 0; n < 3; n++)
        {
          _dR_dF(k, l, m, n) = 0.;
          for (p = 0; p < 3; p++)
            for (q = 0; q < 3; q++)
              _dR_dF(k, l, m, n) +=
                  R(k, p) * (Uhat(p, q) * R(m, q) * Uhat(n, l) - Uhat(p, n) * R(m, q) * Uhat(q, l));

          _dR_dF(k, l, m, n) /= Uhat_det;
        }

  return _dR_dF;
}
