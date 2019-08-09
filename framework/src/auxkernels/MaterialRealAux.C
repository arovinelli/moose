//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MaterialRealAux.h"

registerMooseObject("MooseApp", MaterialRealAux);

template <>
InputParameters
validParams<MaterialRealAux>()
{
  InputParameters params = validParams<MaterialAuxBase<Real>>();
  params.addClassDescription("Outputs element volume-averaged material properties");
  params.addParam<unsigned int>("qp_idx", 0, "output only the selected qp value");
  return params;
}

MaterialRealAux::MaterialRealAux(const InputParameters & parameters)
  : MaterialAuxBase<Real>(parameters),
    _selected_qp_bool(parameters.isParamSetByUser("qp_idx") ? true : false),
    _selected_qp(getParam<unsigned int>("qp_idx"))

{
}

Real
MaterialRealAux::getRealValue()
{
  if (!_selected_qp_bool)
    return _prop[_qp];
  else
    return _prop[_selected_qp];
}
