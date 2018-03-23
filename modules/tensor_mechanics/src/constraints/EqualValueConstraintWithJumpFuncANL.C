//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EqualValueConstraintWithJumpFuncANL.h"
#include "SubProblem.h"
#include "FEProblem.h"
#include "Function.h"

#include <iostream>
using namespace std;

registerMooseObject("TensorMechanicsApp", EqualValueConstraintWithJumpFuncANL);

template <>
InputParameters
validParams<EqualValueConstraintWithJumpFuncANL>()
{
  InputParameters params = validParams<MortarConstraint>();
  params.addRequiredParam<FunctionName>("funcJump", "The forcing function.");
  return params;
}

EqualValueConstraintWithJumpFuncANL::EqualValueConstraintWithJumpFuncANL(const InputParameters & parameters)
  : MortarConstraint(parameters), _funcJump(getFunction("funcJump"))
{
}

Real
EqualValueConstraintWithJumpFuncANL::computeQpResidual()
{
Real jumpValue = _funcJump.value(_t, _q_point[_qp]);
//  cout << "The normal of this quadrature point master";
//  cout << _normals_master[_qp];
//  cout << "The normal of this quadrature point slave";
//  cout << _normals_slave[_qp](0);

  return (_u_master[_qp] - _u_slave[_qp] + jumpValue) * _test[_i][_qp];
}

Real
EqualValueConstraintWithJumpFuncANL::computeQpResidualSide(Moose::ConstraintType res_type)
{
  switch (res_type)
  {
    case Moose::Master:
      return _lambda[_qp] * _test_master[_i][_qp];
    case Moose::Slave:
      return -_lambda[_qp] * _test_slave[_i][_qp];
    default:
      return 0;
  }
}

Real
EqualValueConstraintWithJumpFuncANL::computeQpJacobianSide(Moose::ConstraintJacobianType jac_type)
{
  switch (jac_type)
  {
    case Moose::MasterMaster:
    case Moose::SlaveMaster:
      return _phi[_j][_qp] * _test_master[_i][_qp];

    case Moose::MasterSlave:
    case Moose::SlaveSlave:
      return -_phi[_j][_qp] * _test_slave[_i][_qp];

    default:
      return 0;
  }
}
