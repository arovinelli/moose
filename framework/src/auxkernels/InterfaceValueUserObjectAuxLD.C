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
  params.addClassDescription(
      "Get stored value from the specified UO and save them in an LD element.");

  return params;
}

InterfaceValueUserObjectAuxLD::InterfaceValueUserObjectAuxLD(const InputParameters & parameters)
  : InterfaceValueUserObjectAux(parameters)
{
}

Real
InterfaceValueUserObjectAuxLD::computeValue()
{
  return _interface_uo.getQpValueForLD(_current_elem->id(), _qp);
}
