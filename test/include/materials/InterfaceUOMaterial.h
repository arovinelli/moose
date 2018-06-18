//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef INTERFACEUOMATERIAL_H
#define INTERFACEUOMATERIAL_H

#include "Material.h"
#include "InterfaceUO.h"

class InterfaceUOMaterial;

template <>
InputParameters validParams<InterfaceUOMaterial>();

/**
 *
 */
class InterfaceUOMaterial : public Material
{
public:
  InterfaceUOMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  // virtual void initQpStatefulProperties() override;

  /// User objects computing the displacement jump
  const InterfaceUO & _interface_uo;

  /// the average value of the material property on the interface
  MaterialProperty<Real> & _material_property_average;
  /// the jump in the variable value accros teh interface
  MaterialProperty<Real> & _variable_jump;
};

#endif // INTERFACEUOMATERIAL_H
