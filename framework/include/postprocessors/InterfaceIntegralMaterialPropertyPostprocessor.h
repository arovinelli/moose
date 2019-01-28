//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef INTERFACEINTEGRALMATERIALPROPERTYPOSTPROCESSOR_H
#define INTERFACEINTEGRALMATERIALPROPERTYPOSTPROCESSOR_H

#include "InterfaceIntegralPostprocessor.h"

// Forward Declarations
class InterfaceIntegralMaterialPropertyPostprocessor;

template <>
InputParameters validParams<InterfaceIntegralMaterialPropertyPostprocessor>();

/**
 * This postprocessor computes a volume integral of the specified Material
 * Property.
 *
 * Note that specializations of this integral are possible by deriving from this
 * class and overriding computeQpIntegral().
 */
class InterfaceIntegralMaterialPropertyPostprocessor
    : public InterfaceIntegralPostprocessor {
public:
  InterfaceIntegralMaterialPropertyPostprocessor(
      const InputParameters &parameters);

protected:
  virtual Real computeQpIntegral() override;

  const MaterialProperty<Real> &_mp;
  const MaterialProperty<Real> &_mp_neighbor;
};

#endif
