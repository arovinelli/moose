//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterfaceValueUserObjectAuxLD.h"

registerMooseObject("MooseApp", InterfaceValueUserObjectAuxLD);

template <>
InputParameters
validParams<InterfaceValueUserObjectAuxLD>()
{
  InputParameters params = validParams<InterfaceValueUserObjectAux>();
  params.addRequiredParam<UserObjectName>("LD_map_UO",
                                          "The name of the interface user object to use");
  params.addClassDescription(
      "Get stored value from the specified UO and save them in an LD element.");

  return params;
}

InterfaceValueUserObjectAuxLD::InterfaceValueUserObjectAuxLD(const InputParameters & parameters)
  : InterfaceValueUserObjectAux(parameters), _LDmapUO(getUserObject<map2LDelem>("LD_map_UO"))
{
}

Real
InterfaceValueUserObjectAuxLD::computeValue()
{
  std::pair<dof_id_type, unsigned int> bulk_elem_side = _LDmapUO.getLDNeighbor(_current_elem->id());

  return _interface_uo.getQpValue(bulk_elem_side.first, bulk_elem_side.second, _qp);
}
