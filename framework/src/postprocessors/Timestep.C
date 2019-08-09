//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Timestep.h"
#include "FEProblem.h"

registerMooseObject("MooseApp", Timestep);

template <>
InputParameters
validParams<Timestep>()
{
  InputParameters params = validParams<GeneralPostprocessor>();
  params.addClassDescription("Reports the timestep size");
  return params;
}

Timestep::Timestep(const InputParameters & parameters)
  : GeneralPostprocessor(parameters), _feproblem(dynamic_cast<FEProblemBase &>(_subproblem))
{
}

Real
Timestep::getValue()
{
  return _feproblem.timeStep();
}
