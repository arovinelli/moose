//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterfaceAveragePostprocessor.h"

#include "libmesh/quadrature.h"

template <>
InputParameters
validParams<InterfaceAveragePostprocessor>()
{
  InputParameters params = validParams<InterfaceIntegralPostprocessor>();
  params.addClassDescription(
      "Postprocessor class adding basic capabilities to compute an integral over an interface. "
      "This class is still abstract, refer to InterfaceIntegralVariableValuePostprocessor for a "
      "derived class example");
  return params;
}

InterfaceAveragePostprocessor::InterfaceAveragePostprocessor(const InputParameters & parameters)
  : InterfaceIntegralPostprocessor(parameters)
{
}

Real
InterfaceAveragePostprocessor::getValue()
{
  InterfaceIntegralPostprocessor::getValue();
  return _integral_value / _interface_master_area;
}
