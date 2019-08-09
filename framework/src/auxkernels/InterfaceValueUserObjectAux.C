//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterfaceValueUserObjectAux.h"

registerMooseObject("MooseApp", InterfaceValueUserObjectAux);

template <>
InputParameters
validParams<InterfaceValueUserObjectAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<UserObjectName>("interface_uo_name",
                                          "The name of the interface user object to use");
  params.addClassDescription("Get stored value from the specified InterfaceQpValueUserObjectBase.");
  params.addParam<unsigned int>("qp_idx", 0, "output only the selected qp value");
  return params;
}

InterfaceValueUserObjectAux::InterfaceValueUserObjectAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _interface_uo(getUserObject<InterfaceQpValueUserObjectBase>("interface_uo_name")),
    _selected_qp_bool(parameters.isParamSetByUser("qp_idx") ? true : false),
    _selected_qp(getParam<unsigned int>("qp_idx"))
{
}

Real
InterfaceValueUserObjectAux::computeValue()
{
  unsigned int my_qp = _qp;
  if (_selected_qp_bool)
    my_qp = _selected_qp;

  return _interface_uo.getQpValue(_current_elem->id(), _current_side, my_qp);
}
