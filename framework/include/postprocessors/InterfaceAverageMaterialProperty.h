//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef INTERFACEAVERAGEMATERIALPROPERTY_H
#define INTERFACEAVERAGEMATERIALPROPERTY_H

#include "InterfaceIntegralMaterialPropertyPostprocessor.h"

// Forward Declarations
class InterfaceAverageMaterialProperty;

template <> InputParameters validParams<InterfaceAverageMaterialProperty>();

/**
 * This postprocessor computes a volume integral of the specified variable.
 *
 * Note that specializations of this integral are possible by deriving from this
 * class and overriding computeQpIntegral().
 */
class InterfaceAverageMaterialProperty
    : public InterfaceIntegralMaterialPropertyPostprocessor {
public:
  InterfaceAverageMaterialProperty(const InputParameters &parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual Real getValue() override;
  virtual void threadJoin(const UserObject &y) override;

protected:
  virtual Real volume();
  Real _volume;
};

#endif
