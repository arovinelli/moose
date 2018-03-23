//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef EQUALVALUECONSTRAINTWITHCOUPLEDVAR_H
#define EQUALVALUECONSTRAINTWITHCOUPLEDVAR_H

#include "MortarConstraint.h"

class EqualValueConstraintWithCoupledVarANL;

template <>
InputParameters validParams<EqualValueConstraintWithCoupledVarANL>();

/**
 * Constrain the value of a variable to be the same on both sides of an
 * interface.
 */
class EqualValueConstraintWithCoupledVarANL : public MortarConstraint
{
public:
  EqualValueConstraintWithCoupledVarANL(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpResidualSide(Moose::ConstraintType res_type) override;
  virtual Real computeQpJacobianSide(Moose::ConstraintJacobianType jac_type) override;


  unsigned int _coupled_var_ID; // coupled variable identifier
  Function & _funcJump;
};

#endif /* EQUALVALUECONSTRAINTWITHCOUPLEDVARANL_H */
