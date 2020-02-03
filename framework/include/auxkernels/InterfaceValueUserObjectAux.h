//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"
#include "InterfaceQpUserObjectBase.h"

// Forward Declarations
class InterfaceValueUserObjectAux;

template <>
InputParameters validParams<InterfaceValueUserObjectAux>();

/**
 * AuxKernel creating an AuxVariable from values stored in an InterfaceQpUserObjectBase
 */
class InterfaceValueUserObjectAux : public AuxKernel
{
public:
  static InputParameters validParams();

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  InterfaceValueUserObjectAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
  const InterfaceQpUserObjectBase & _interface_uo;
  const bool _return_side_average;
};
