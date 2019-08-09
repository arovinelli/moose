//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef INTERFACEUSEROBJECTQPAUXLD_H
#define INTERFACEUSEROBJECTQPAUXLD_H

#include "InterfaceQpValueUserObject.h"
#include "map2LDelem.h"

// Forward Declarations
class InterfaceValueUserObjectAuxLD;

template <>
InputParameters validParams<InterfaceValueUserObjectAuxLD>();

/**
 * AuxKernel creating an AuxVariable from values stored in an InterfaceQpValueUserObject
 */
class InterfaceValueUserObjectAuxLD : public InterfaceValueUserObjectAux
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  InterfaceValueUserObjectAuxLD(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
  const map2LDelem & _LDmapUO;
};

#endif // INTERFACEUSEROBJECTQPAUXLD_H
