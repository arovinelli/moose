//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GenericConstantMaterial.h"
#include "petscsnes.h"
// Forward Declarations
class Test_PETSC;
template <>
InputParameters validParams<Test_PETSC>();

/**
 * This material automatically declares as material properties whatever is
 * passed to it through the parameters 'prop_names' and uses the values from
 * 'prop_values' as the values for those properties.
 *
 * This is not meant to be used in a production capacity... and instead is meant
 * to be used during development phases for ultimate flexibility.
 */
class Test_PETSC : public GenericConstantMaterial
{
public:
  Test_PETSC(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;
  // const Real _f = 50;
  SNES _petsc_snes;       /* nonlinear solver context */
  KSP _petsc_ksp;         /* linear solver context */
  PC _petsc_pc;           /* preconditioner context */
  Vec _petsc_x, _petsc_r; /* solution, residual vectors */
  Mat _petsc_J;           /* Jacobian matrix */
  PetscErrorCode _petsc_ierr;
  PetscErrorCode * _petsc_ierr_pt;
  PetscScalar _pfive = .5, *xx;
};
