//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearNonLinearIterationMaterial.h"

registerMooseObject("MooseTestApp", LinearNonLinearIterationMaterial);

InputParameters
LinearNonLinearIterationMaterial::validParams()
{
  InputParameters params = GenericConstantMaterial::validParams();
  params.addClassDescription("Material whose property  is equal to t_step+n_linear*n_nonlinear");
  return params;
}

LinearNonLinearIterationMaterial::LinearNonLinearIterationMaterial(
    const InputParameters & parameters)
  : GenericConstantMaterial(parameters), _mat_prop(declareProperty<Real>("mat_prop"))
{
}

void
LinearNonLinearIterationMaterial::computeQpProperties()
{
  _mat_prop[_qp] = _t_step + _fe_problem.nNonlinearIterations() * _fe_problem.nLinearIterations();
}
