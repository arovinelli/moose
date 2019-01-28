//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterfaceAverageValue.h"

registerMooseObject("MooseApp", InterfaceAverageValue);

template <>
InputParameters
validParams<InterfaceAverageValue>()
{
  InputParameters params = validParams<InterfaceIntegralVariablePostprocessor>();
  params.addClassDescription("Computes the average value of a variable on a "
                             "sideset. Note that this cannot be used on the "
                             "centerline of an axisymmetric model.");
  return params;
}

InterfaceAverageValue::InterfaceAverageValue(const InputParameters & parameters)
  : InterfaceIntegralVariablePostprocessor(parameters), _volume(0)
{
}

void
InterfaceAverageValue::initialize()
{
  InterfaceIntegralVariablePostprocessor::initialize();
  _volume = 0;
}

void
InterfaceAverageValue::execute()
{
  InterfaceIntegralVariablePostprocessor::execute();
  _volume += volume();
}

Real
InterfaceAverageValue::getValue()
{
  Real integral = InterfaceIntegralVariablePostprocessor::getValue();
  gatherSum(_volume);
  return integral / _volume;
}

Real
InterfaceAverageValue::volume()
{
  return _current_side_volume;
}

void
InterfaceAverageValue::threadJoin(const UserObject & y)
{
  InterfaceIntegralVariablePostprocessor::threadJoin(y);
  const InterfaceAverageValue & pps = static_cast<const InterfaceAverageValue &>(y);
  _volume += pps._volume;
}
