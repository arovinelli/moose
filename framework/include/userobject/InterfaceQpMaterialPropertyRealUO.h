//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceQpMaterialPropertyUserObjectBase.h"

// Forward declarations
class InterfaceQpMaterialPropertyRealUO;

template <>
InputParameters validParams<InterfaceQpMaterialPropertyRealUO>();

/**
 * This userobject collect values of a variable across an interface for each QP and compute a
 * scalar. The computed scalar value depends on the given parameter _interface_value_type\
 * _interface_value_type (see IntervafeValueTools).
 */
class InterfaceQpMaterialPropertyRealUO : public InterfaceQpMaterialPropertyUserObjectBase<Real>
{

public:
  static InputParameters validParams();
  /**
   * Class constructor
   * @param parameters The input parameters for this object
   */
  InterfaceQpMaterialPropertyRealUO(const InputParameters & parameters);
  virtual ~InterfaceQpMaterialPropertyRealUO(){};

protected:
  virtual Real computeValueMaster(const unsigned int /*qp*/) override;
  virtual Real computeValueSlave(const unsigned int /*qp*/) override;
};
