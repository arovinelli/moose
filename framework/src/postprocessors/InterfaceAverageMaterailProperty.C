//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterfaceAverageMaterialProperty.h"

registerMooseObject("MooseApp", InterfaceAverageMaterialProperty);

template <>
InputParameters
validParams<InterfaceAverageMaterialProperty>()
{
  InputParameters params = validParams<InterfaceIntegralMaterialPropertyPostprocessor>();
  params.addClassDescription("Computes the average value of a material property on a "
                             "interface. Note that this cannot be used on the "
                             "centerline of an axisymmetric model.");
  return params;
}

InterfaceAverageMaterialProperty::InterfaceAverageMaterialProperty(
    const InputParameters & parameters)
  : InterfaceIntegralMaterialPropertyPostprocessor(parameters), _volume(0)
{
}

void
InterfaceAverageMaterialProperty::initialize()
{
  InterfaceIntegralMaterialPropertyPostprocessor::initialize();
  _volume = 0;
}

void
InterfaceAverageMaterialProperty::execute()
{
  InterfaceIntegralMaterialPropertyPostprocessor::execute();
  _volume += volume();
}

Real
InterfaceAverageMaterialProperty::getValue()
{
  Real integral = InterfaceIntegralMaterialPropertyPostprocessor::getValue();
  gatherSum(_volume);
  return integral / _volume;
}

Real
InterfaceAverageMaterialProperty::volume()
{
  return _current_side_volume;
}

void
InterfaceAverageMaterialProperty::threadJoin(const UserObject & y)
{
  InterfaceIntegralMaterialPropertyPostprocessor::threadJoin(y);
  const InterfaceAverageMaterialProperty & pps =
      static_cast<const InterfaceAverageMaterialProperty &>(y);
  _volume += pps._volume;
}
