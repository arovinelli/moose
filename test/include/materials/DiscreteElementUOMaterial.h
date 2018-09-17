//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef DISCRETEELEMENTUOMATERIAL
#define DISCRETEELEMENTUOMATERIAL

#include "Material.h"
#include "DiscreteElementWritingMP.h"

// Forward Declarations
class DiscreteElementUOMaterial;

template <>
InputParameters validParams<DiscreteElementUOMaterial>();

/**
 * Stateful material class that defines a few properties.
 */
class DiscreteElementUOMaterial : public Material
{
public:
  DiscreteElementUOMaterial(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();

private:
  /// User objects that writing/returning updated MP
  const DiscreteElementWritingMP & _discrete_elem_uo;
  const unsigned int _n_uo_mp;
  // MP instantiated by the UOs
  std::vector<MaterialProperty<std::vector<Real>> *> _uomp_values;
  std::vector<const MaterialProperty<std::vector<Real>> *> _uomp_values_old;

  // MaterialProperty<Real> & _test_MP_1;
  // const MaterialProperty<Real> & _test_MP_1_old;
  // MaterialProperty<Real> & _test_MP_2;
  // const MaterialProperty<Real> & _test_MP_2_old;
};

#endif // DISCRETEELEMENTUOMATERIAL
