//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef DISCRETEELEMENTWRITINGMP_H
#define DISCRETEELEMENTWRITINGMP_H

#include "DiscreteElementUserObject.h"

class DiscreteElementWritingMP;

template <>
InputParameters validParams<DiscreteElementWritingMP>();

/**
 * Crystal plasticity system userobject base class.
 */
class DiscreteElementWritingMP : public DiscreteElementUserObject
{
public:
  DiscreteElementWritingMP(const InputParameters & parameters);

  /// Returns Material Properties Sizes
  virtual unsigned int getMaterialPropertySize(unsigned int /*mp_index*/) const;

  /// Returns the Material Properties Names
  virtual std::string getMaterialPropertyName(unsigned int /*mp_index*/) const;

  /// Returns the statefuln character of all MP
  virtual bool getMaterialPropertyStateful(unsigned int /*mp_index*/) const;

  /// Returns the Material Properties Names
  virtual std::vector<Real> getInitialMaterialPropertyValues(unsigned int /*mp_index*/) const;

  /// Returns the Material Properties Names
  virtual unsigned int getNumberMaterialProperty() const;

protected:
  const std::vector<std::string> _MP_names;
  const std::vector<unsigned int> _MP_sizes;
  const std::vector<unsigned int> _MP_stateful;
  const std::vector<std::vector<Real>> _MP_initial_values;
  const unsigned int _n_MP;
  std::vector<std::vector<Real>> ResizeInitialValues() const;
};

#endif // DISCRETEELEMENTWRITINGMP_H
