//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterfaceIntegralPostprocessor.h"

#include "libmesh/quadrature.h"

template <> InputParameters validParams<InterfaceIntegralPostprocessor>() {
  InputParameters params = validParams<InterfacePostprocessor>();
  params.addParam<MooseEnum>(
      "integral_type", InterfaceIntegralPostprocessor::integralTypeOptions(),
      "the type fo integral we want calcualte. Options "
      "are: average, jump_master_slave, "
      "jump_slave_master, jump_abs, master, slave");
  return params;
}

InterfaceIntegralPostprocessor::InterfaceIntegralPostprocessor(
    const InputParameters &parameters)
    : InterfacePostprocessor(parameters), _qp(0), _integral_value(0),
      _integral_type(parameters.get<MooseEnum>("integral_type")) {}

void InterfaceIntegralPostprocessor::initialize() { _integral_value = 0; }

void InterfaceIntegralPostprocessor::execute() {
  _integral_value += computeIntegral();
}

Real InterfaceIntegralPostprocessor::getValue() {
  gatherSum(_integral_value);
  return _integral_value;
}

void InterfaceIntegralPostprocessor::threadJoin(const UserObject &y) {
  const InterfaceIntegralPostprocessor &pps =
      static_cast<const InterfaceIntegralPostprocessor &>(y);
  _integral_value += pps._integral_value;
}

Real InterfaceIntegralPostprocessor::computeIntegral() {
  Real sum = 0;
  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    sum += _JxW[_qp] * _coord[_qp] * computeQpIntegral();
  return sum;
}

Real InterfaceIntegralPostprocessor::computeIntegralType(
    const Real value_master, const Real value_slave) {
  Real result;

  switch (_integral_type) {
  case 0: /*average*/
    result = (value_master + value_slave) * 0.5;
    break;
  case 1: /*jump_master_minus_slave*/
    result = (value_master - value_slave);
    break;
  case 2: /*jump_slave_minus_master*/
    result = (value_slave - value_master);
    break;
  case 3: /*jump_abs*/
    result = std::abs(value_slave - value_master);
    break;
  case 4: /*master*/
    result = value_master;
    break;
  case 5: /*slave*/
    result = value_slave;
    break;
  default:
    mooseError(
        "InterfaceIntegralMaterialPropertyPostprocessor: the supplied intergal "
        "type is not in the list. Avaialble options are: ",
        _integral_type.getRawNames());
  }
  return result;
}

MooseEnum InterfaceIntegralPostprocessor::integralTypeOptions() {
  return MooseEnum("average jump_master_minus_slave jump_slave_minus_master "
                   "jump_abs master slave",
                   "average");
}
