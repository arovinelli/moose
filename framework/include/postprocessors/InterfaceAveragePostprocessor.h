//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceIntegralPostprocessor.h"

// Forward Declarations
class InterfaceAveragePostprocessor;

template <>
InputParameters validParams<InterfaceAveragePostprocessor>();

/**
 * This postprocessor add generel capabilities to the InterfaceIntegralPostprocessor to compute an
 * integral over an interface. To actually compute an integral one must derive from this class,
 * specialize computeQpIntegral() and give access to either a varaible or a material property
 **/
class InterfaceAveragePostprocessor : public InterfaceIntegralPostprocessor
{
public:
  InterfaceAveragePostprocessor(const InputParameters & parameters);

protected:
  Real getValue() override;
  Real computeQpIntegral() override = 0; // MUST BE OVERRIDDEN
};
