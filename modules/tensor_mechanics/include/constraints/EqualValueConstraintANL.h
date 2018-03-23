//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef EQUALVALUECONSTRAINTANL_H
#define EQUALVALUECONSTRAINTANL_H

#include "MortarConstraint.h"

class EqualValueConstraintANL;

template <>
InputParameters validParams<EqualValueConstraintANL>();

/**
 * Constrain the value of a variable to be the same on both sides of an
 * interface.
 */
class EqualValueConstraintANL : public MortarConstraint
{
public:
  EqualValueConstraintANL(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpResidualSide(Moose::ConstraintType res_type) override;
  virtual Real computeQpJacobianSide(Moose::ConstraintJacobianType jac_type) override;
};

#endif /* EQUALVALUECONSTRAINTANL_H */
