//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DiscreteElementUOMaterial.h"

registerMooseObject("MooseTestApp", DiscreteElementUOMaterial);

template <>
InputParameters
validParams<DiscreteElementUOMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<UserObjectName>("discrete_elem_uo",
                                          "the name of the discrete element user object");
  return params;
}

DiscreteElementUOMaterial::DiscreteElementUOMaterial(const InputParameters & parameters)
  : Material(parameters),

    _discrete_elem_uo(getUserObjectByName<DiscreteElementWritingMP>(
        parameters.get<UserObjectName>("discrete_elem_uo"))),
    _n_uo_mp(_discrete_elem_uo.getNumberMaterialProperty())

{
  // initialize the userobject material property container
  _uomp_values.resize(_n_uo_mp);
  _uomp_values_old.resize(_n_uo_mp);

  for (unsigned int mp_index = 0; mp_index < _n_uo_mp; mp_index++)
  {
    // declare a material property
    _uomp_values[mp_index] =
        &declareProperty<std::vector<Real>>(_discrete_elem_uo.getMaterialPropertyName(mp_index));

    // for a stateful material property get the old value
    if (_discrete_elem_uo.getMaterialPropertyStateful(mp_index))
      _uomp_values_old[mp_index] = &getMaterialPropertyOld<std::vector<Real>>(
          _discrete_elem_uo.getMaterialPropertyName(mp_index));
  }
}

void
DiscreteElementUOMaterial::initQpStatefulProperties()
{
  for (unsigned int mp_index = 0; mp_index < _n_uo_mp; mp_index++)
  {

    // check if the i-th material property is stateful
    if (_discrete_elem_uo.getMaterialPropertyStateful(mp_index))
    {
      // resize the material property
      (*_uomp_values[mp_index])[_qp].resize(_discrete_elem_uo.getMaterialPropertySize(mp_index));

      // assign intial values to the stateful material property
      for (unsigned int i = 0; i < _discrete_elem_uo.getMaterialPropertySize(mp_index); i++)
      {
        (*_uomp_values[mp_index])[_qp] =
            _discrete_elem_uo.getInitialMaterialPropertyValues(mp_index);
      }
    }
  }
}

void
DiscreteElementUOMaterial::computeQpProperties()
{

  for (unsigned int mp_index = 0; mp_index < _n_uo_mp; mp_index++)
  {
    // check if the i-th material property is stateful
    if (!_discrete_elem_uo.getMaterialPropertyStateful(mp_index))
      // if not stateful resize
      (*_uomp_values[mp_index])[_qp].resize(_discrete_elem_uo.getMaterialPropertySize(mp_index));

    for (unsigned int i = 0; i < _discrete_elem_uo.getMaterialPropertySize(mp_index); i++)
    {
      if (_discrete_elem_uo.getMaterialPropertyStateful(mp_index))
      {
        (*_uomp_values[mp_index])[_qp][i] =
            _discrete_elem_uo.getInitialMaterialPropertyValues(mp_index)[i] * _dt +
            (*_uomp_values_old[mp_index])[_qp][i];
        if (_current_elem->id() == 0 && _qp == 0)
          std::cout << "dt " << _dt << " MP_name "
                    << _discrete_elem_uo.getMaterialPropertyName(mp_index) << " MP_val "
                    << (*_uomp_values[mp_index])[_qp][i] << " MP_val_old "
                    << (*_uomp_values_old[mp_index])[_qp][i] << std::endl;
      }

      else
        (*_uomp_values[mp_index])[_qp][i] =
            _discrete_elem_uo.getInitialMaterialPropertyValues(mp_index)[i] * 2.;
    }
  }
}
