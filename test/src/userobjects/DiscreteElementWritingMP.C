//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DiscreteElementWritingMP.h"

registerMooseObject("MooseTestApp", DiscreteElementWritingMP);

template <>
InputParameters
validParams<DiscreteElementWritingMP>()
{
  InputParameters params = validParams<DiscreteElementUserObject>();
  params.addRequiredParam<std::vector<std::string>>("MP_names", "The MPs' names");
  params.addRequiredParam<std::vector<unsigned int>>("MP_sizes",
                                                     "A vector containing the MPs' sizes");
  params.addRequiredParam<std::vector<unsigned int>>(
      "MP_stateful", "flag specifying if the i-th property is stateful");
  params.addRequiredParam<std::vector<Real>>(
      "MP_initial_values", "The MP's intialia values, all in a singel vector and ordered by name");
  params.addClassDescription("Discrete Element UO Writing the updated value of a MP");

  return params;
}

DiscreteElementWritingMP::DiscreteElementWritingMP(const InputParameters & parameters)
  : DiscreteElementUserObject(parameters),
    _MP_names(getParam<std::vector<std::string>>("MP_names")),
    _MP_sizes(getParam<std::vector<unsigned int>>("MP_sizes")),
    _MP_stateful(getParam<std::vector<unsigned int>>("MP_stateful")),
    _MP_initial_values(ResizeInitialValues()),
    _n_MP(_MP_names.size())
{
}

unsigned int
DiscreteElementWritingMP::getNumberMaterialProperty() const
{
  return _n_MP;
}

std::string
DiscreteElementWritingMP::getMaterialPropertyName(unsigned int mp_index) const
{
  return _MP_names[mp_index];
}

unsigned int
DiscreteElementWritingMP::getMaterialPropertySize(unsigned int mp_index) const
{
  return _MP_sizes[mp_index];
}

bool
DiscreteElementWritingMP::getMaterialPropertyStateful(unsigned int mp_index) const
{
  return _MP_stateful[mp_index];
}

std::vector<Real>
DiscreteElementWritingMP::getInitialMaterialPropertyValues(unsigned int mp_index) const
{
  return _MP_initial_values[mp_index];
}

std::vector<std::vector<Real>>
DiscreteElementWritingMP::ResizeInitialValues() const
{
  std::vector<std::vector<Real>> temp;
  std::vector<Real> temp_init_values = getParam<std::vector<Real>>("MP_initial_values");
  unsigned int n_subvector = _MP_sizes.size();
  temp.resize(n_subvector);
  unsigned int c = 0;
  for (unsigned int i = 0; i < n_subvector; i++)
  {
    unsigned int subvector_size = _MP_sizes[i];
    temp[i].resize(subvector_size);
    for (unsigned int j = 0; j < subvector_size; j++)
    {
      temp[i][j] = temp_init_values[c];
      c += 1;
    }
  }
  return temp;
}
