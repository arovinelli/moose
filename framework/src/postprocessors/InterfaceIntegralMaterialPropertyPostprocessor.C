//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterfaceIntegralMaterialPropertyPostprocessor.h"

registerMooseObject("MooseApp", InterfaceIntegralMaterialPropertyPostprocessor);

template <>
InputParameters
validParams<InterfaceIntegralMaterialPropertyPostprocessor>()
{
  InputParameters params = validParams<InterfaceIntegralPostprocessor>();

  params.addRequiredParam<MaterialPropertyName>(
      "property", "The name of the material property on the master side of the interface");
  params.addParam<MaterialPropertyName>(
      "neighbor_property",
      "The name of the material property on the slave side of the interface.  "
      "Default to the master material property name if omitted");
  params.addParam<bool>("use_old_prop",
                        false,
                        "whetever to use the current or previous timestep "
                        "material property value");
  return params;
}

InterfaceIntegralMaterialPropertyPostprocessor::InterfaceIntegralMaterialPropertyPostprocessor(
    const InputParameters & parameters)
  : InterfaceIntegralPostprocessor(parameters),
    _mp(parameters.get<bool>("use_old_prop") ? getMaterialPropertyOld<Real>("property")
                                             : getMaterialProperty<Real>("property")),
    _mp_neighbor(parameters.isParamSetByUser("neighbor_property")
                     ? (parameters.get<bool>("use_old_prop")
                            ? getNeighborMaterialPropertyOld<Real>("neighbor_property")
                            : getNeighborMaterialProperty<Real>("neighbor_property"))
                     : (parameters.get<bool>("use_old_prop")
                            ? getNeighborMaterialPropertyOld<Real>("property")
                            : getNeighborMaterialProperty<Real>("property")))
{
}

Real
InterfaceIntegralMaterialPropertyPostprocessor::computeQpIntegral()
{
  return computeIntegralType(_mp[_qp], _mp_neighbor[_qp]);
}
