//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MaterialRealVectorValueAux.h"

registerMooseObject("MooseApp", MaterialRealVectorValueAux);

template <>
InputParameters
validParams<MaterialRealVectorValueAux>()
{
  InputParameters params = validParams<MaterialAuxBase<>>();
  params.addParam<unsigned int>("component", 0, "The vector component to consider for this kernel");
  params.addParam<unsigned int>("qp_idx", 0, "output only the selected qp value");

  return params;
}

MaterialRealVectorValueAux::MaterialRealVectorValueAux(const InputParameters & parameters)
  : MaterialAuxBase<RealVectorValue>(parameters),
    _component(getParam<unsigned int>("component")),
    _selected_qp_bool(parameters.isParamSetByUser("qp_idx") ? true : false),
    _selected_qp(getParam<unsigned int>("qp_idx"))
{
  if (_component > LIBMESH_DIM)
    mooseError(
        "The component ", _component, " does not exist for ", LIBMESH_DIM, " dimensional problems");
}

Real
MaterialRealVectorValueAux::getRealValue()
{
  if (!_selected_qp_bool)
    return _prop[_qp](_component);
  else
    return _prop[_selected_qp](_component);
}
