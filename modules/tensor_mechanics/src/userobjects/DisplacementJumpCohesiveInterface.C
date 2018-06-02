//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DisplacementJumpCohesiveInterface.h"
#include "MooseError.h"
// #include "RankTwoTensor.h"
#include "RotationMatrix.h"

registerMooseObject("TensorMechanicsApp", DisplacementJumpCohesiveInterface);

template <>
InputParameters
validParams<DisplacementJumpCohesiveInterface>()
{
  InputParameters params = validParams<InternalSideUserObject>();

  params.addClassDescription("Compute the dispalcment jump across a cohesive interface");
  params.addRequiredCoupledVar("disp_x",
                               "variable containing the X component of"
                               "the displacement on the master side");
  params.addRequiredCoupledVar("disp_x_neighbor",
                               "variable containing the X component"
                               " of the displacement on the slave side");
  params.addCoupledVar("disp_y",
                       "variable containing the Y component of"
                       "the displacement on the master side");
  params.addCoupledVar("disp_y_neighbor",
                       "variable containing the Y "
                       "component of the displacement on the slave side");
  params.addCoupledVar("disp_z",
                       "variable containing the Z component of"
                       "the displacement on the master side");
  params.addCoupledVar("disp_z_neighbor",
                       "variable containing the Z "
                       "component of the displacement on the slave side");

  return params;
}

DisplacementJumpCohesiveInterface::DisplacementJumpCohesiveInterface(const InputParameters & params)
  : InternalSideUserObject(params),
    _disp_x(coupledValue("disp_x")),
    _disp_x_neighbor(coupledNeighborValue("disp_x_neighbor")),
    _disp_y(_mesh.dimension() >= 2 ? coupledValue("disp_y") : _zero),
    _disp_y_neighbor(_mesh.dimension() >= 2 ? coupledNeighborValue("disp_y_neighbor") : _zero),
    _disp_z(_mesh.dimension() >= 3 ? coupledValue("disp_z") : _zero),
    _disp_z_neighbor(_mesh.dimension() >= 3 ? coupledNeighborValue("disp_z_neighbor") : _zero)
{
}

void
DisplacementJumpCohesiveInterface::computeResidaulAndJacobianCoefficients(
    dof_id_type elem,
    unsigned int side,
    unsigned int qp,
    RealVectorValue & Traction,
    RankTwoTensor & TractionSpatialDerivative) const
{
  RealVectorValue Jump;
  RealVectorValue side_normal;
  RealVectorValue JumpLocal;
  RealVectorValue TractionLocal;
  RealTensorValue RotationGlobal2Local;
  RankTwoTensor TractionSpatialDerivativeLocal;

  // std::cout << "RetrieveDisplacementJump UO" << std::endl;
  auto data = _map_values.find(std::make_pair(elem, side));
  if (data != _map_values.end())
  {
    Jump = data->second[qp][0];
    side_normal = data->second[qp][1];

    moveToLocalFrame(Jump, JumpLocal, side_normal, RotationGlobal2Local);

    computeTractionLocal(JumpLocal, TractionLocal);

    computeTractionSpatialDerivativeLocal(JumpLocal, TractionSpatialDerivativeLocal);

    moveBackToGlobalFrame(TractionLocal,
                          TractionSpatialDerivativeLocal,
                          Traction,
                          TractionSpatialDerivative,
                          RotationGlobal2Local);
  }
  else
  {
    mooseError("can't find the given qp");
  }
}

RealVectorValue
DisplacementJumpCohesiveInterface::computeDisplacementJump(unsigned int qp)
{
  RealVectorValue Jump;

  Jump(0) = _disp_x_neighbor[qp] - _disp_x[qp];
  Jump(1) = _disp_y_neighbor[qp] - _disp_y[qp];
  Jump(2) = _disp_z_neighbor[qp] - _disp_z[qp];

  return Jump;
}

// ovverride standard UO functions
void
DisplacementJumpCohesiveInterface::initialize()
{

  // std::cout << "init UO" << std::endl;

  std::vector<dof_id_type> el;
  std::vector<unsigned short int> sl;
  std::vector<boundary_id_type> il;

  std::set<BoundaryID> boundaryList;
  boundaryList.insert(100);
  /*boundaryIDs()*/;

  _mesh.buildSideList(el, sl, il);
  _map_values.clear();

  // std::cout << "boundary list created" << std::endl;
  for (unsigned int i = 0; i < il.size(); i++)
  {
    if (boundaryList.find(il[i]) != boundaryList.end())
    {
      std::pair<dof_id_type, unsigned int> elem_side_pair = std::make_pair(el[i], sl[i]);
      std::vector<std::vector<RealVectorValue>> var_values(0, std::vector<RealVectorValue>(2));
      _map_values[elem_side_pair] = var_values;
    }
  }
  // std::cout << "init UO completed" << std::endl;
}

void
DisplacementJumpCohesiveInterface::execute()
{

  auto it = _map_values.find(std::make_pair(_current_elem->id(), _current_side));
  if (it != _map_values.end())
  {
    // std::cout << " EL_ID " << _current_elem->id() << std::endl;
    // std::cout << "    SIDE " << _current_side << std::endl;
    auto & vec = _map_values[std::make_pair(_current_elem->id(), _current_side)];
    vec.resize(_qrule->n_points());
    for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
    {
      vec[qp].resize(2, RealVectorValue());
      // std::cout << "         qp " << qp << std::endl;

      vec[qp][0] = computeDisplacementJump(qp);
      vec[qp][1] = _normals[qp];

      // std::cout << "              Jump: " << vec[qp][0](0) << " " << vec[qp][0](1) << " "
      //           << vec[qp][0](2) << std::endl;
      // std::cout << "              Normal: " << vec[qp][1](0) << " " << vec[qp][1](1) << " "
      //           << vec[qp][1](2) << std::endl;
    }
  }
}

void
DisplacementJumpCohesiveInterface::finalize()
{
}

void
DisplacementJumpCohesiveInterface::threadJoin(const UserObject &)
{
}

void
DisplacementJumpCohesiveInterface::moveToLocalFrame(RealVectorValue & Jump,
                                                    RealVectorValue & JumpLocal,
                                                    RealVectorValue & normals,
                                                    RealTensorValue & RotationGlobal2Local) const
{
  // this is a rotation matrix that will rotate _n to the "x" axis such that the
  // the first compoenent in the local frame represent the opening displacment
  RotationGlobal2Local = RotationMatrix::rotVec1ToVec2(normals, RealVectorValue(1, 0, 0));

  // compute the jump in the lcoal coordinate system
  for (unsigned int i = 0; i < 3; i++)
  {
    JumpLocal(i) = 0; // just to be sure that it reset at each iteration
    for (unsigned int j = 0; j < 3; j++)
    {
      JumpLocal(i) += RotationGlobal2Local(i, j) * Jump(j);
    }
  }
}

void
DisplacementJumpCohesiveInterface::moveBackToGlobalFrame(
    RealVectorValue & TractionLocal,
    RankTwoTensor & TractionSpatialDerivativeLocal,
    RealVectorValue & Traction,
    RankTwoTensor & TractionSpatialDerivative,
    RealTensorValue & RotationGlobal2Local) const
{

  RealTensorValue RotationLocal2Global = RotationGlobal2Local.transpose();
  // rotate traction in the global frame
  for (unsigned int i = 0; i < 3; i++)
  {
    Traction(i) = 0;
    for (unsigned int j = 0; j < 3; j++)
    {
      Traction(i) += RotationLocal2Global(i, j) * TractionLocal(j);
    }
  }

  // rotate traction derivatives in the global frame
  TractionSpatialDerivative = TractionSpatialDerivativeLocal;
  TractionSpatialDerivative.rotate(RotationLocal2Global);
}

void
DisplacementJumpCohesiveInterface::computeTractionLocal(RealVectorValue & _JumpLocal,
                                                        RealVectorValue & TractionLocal) const
{
  std::vector<Real> _deltaU0(3, 0);
  _deltaU0[0] = 1;
  _deltaU0[1] = 0.5;
  _deltaU0[2] = 0.5;

  std::vector<Real> _maxAllowableTraction(3, 0);
  _maxAllowableTraction[0] = 100;
  _maxAllowableTraction[1] = 50;
  _maxAllowableTraction[2] = 50;

  // convention N, T, S
  Real temp, X, expX, A_i, B_i;

  X = 0;
  for (unsigned int i = 0; i < 3; i++)
  {
    temp = _JumpLocal(i) / _deltaU0[i];
    if (i > 0)
    {
      temp *= temp; // square for shear component
    };
    X += temp;
  }

  expX = std::exp(-X);

  for (unsigned int i = 0; i < 3; i++)
  {

    if (i == 0)
    {
      temp = std::exp(1);
    }
    else
    {
      temp = std::sqrt(2 * std::exp(1));
    }
    A_i = _maxAllowableTraction[i] * temp;

    B_i = _JumpLocal(i) / _deltaU0[i];

    TractionLocal(i) = A_i * B_i * expX;
  }

  return;
}

void
DisplacementJumpCohesiveInterface::computeTractionSpatialDerivativeLocal(
    RealVectorValue & _JumpLocal, RankTwoTensor & TractionDerivativeLocal) const
{

  std::vector<Real> _deltaU0(3, 0);
  _deltaU0[0] = 1;
  _deltaU0[1] = 0.5;
  _deltaU0[2] = 0.5;

  std::vector<Real> _maxAllowableTraction(3, 0);
  _maxAllowableTraction[0] = 100;
  _maxAllowableTraction[1] = 50;
  _maxAllowableTraction[2] = 50;

  // this function compute partial derivates of Tn[0][:], Tt[1][:], Ts[2][:]
  // w.r.t. dun, dut, dus
  // T_i = A_i*B_i*exp(-X) with:
  // A_i = \sigma_i,max * (\alpha_i*e)^{1/\alpha_i} with \alpha_i = 1 for i==n
  // \alpha_i = 2 for i!=n
  // B_i = \delta_u,i / \delta_0,i
  // X = sum_i=1^3{(\delta_u,i / \delta_0,i)^\alpha_i}  with \alpha_i = 1 for i==n
  // \alpha_i = 2 for i!=n

  // dTi_duj = A_i * ( dBi_duj * exp(-X) + B_i * exp(-X) * dX_duj  )
  //         = A_i * ( exp(-X) * (dBi_duj + B_i * dX_duj ) )

  // convention N, T, S
  unsigned int i, j;
  Real expX, temp, X;

  // compute X and the exponential term
  temp = 0;
  X = 0;
  for (i = 0; i < 3; i++)
  {
    temp = _JumpLocal(i) / _deltaU0[i];
    if (i > 0)
      temp *= temp;
    X += temp;
  }
  expX = std::exp(-X);

  // compute partial derivatives in local coordaintes w.r.t. the master surface siplacement
  //            | dTn/dun dTn/dut dTn/dus |
  // dTi_duj  = | dTt/dun dTt/dut dTt/dus | = _TractionDerivativeLocal[i][j]
  //            | dTs/dun dTs/dut dTs/dus |
  Real A_i, B_i;
  Real dBi_dui, dX_duj;

  for (i = 0; i < 3; i++)
  {

    // compute A_i
    if (i == 0) // alpha = 1
      A_i = std::exp(1);
    else // alpha = 2
      A_i = std::sqrt(2 * std::exp(1));

    A_i *= _maxAllowableTraction[i];

    // compute B_i
    B_i = _JumpLocal(i) / _deltaU0[i];

    for (j = 0; j < 3; j++)
    {

      // add term for diagonal entry dBi_dui
      dBi_dui = 0;
      if (i == j)
      {
        dBi_dui = 1 / _deltaU0[j];
      }

      // compute the derivative of the argument of exponential
      if (j == 0) // alpha = 1
        dX_duj = 1. / _deltaU0[j];
      else // alpha = 2
        dX_duj = 2. * _JumpLocal(j) / (_deltaU0[j] * _deltaU0[j]);

      TractionDerivativeLocal(i, j) =
          A_i * expX * (dBi_dui - B_i * dX_duj); // the minus sign is due to exp(-X)
    }
  }

  return;
}
