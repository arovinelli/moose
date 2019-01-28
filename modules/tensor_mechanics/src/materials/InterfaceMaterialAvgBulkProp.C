//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterfaceMaterialAvgBulkProp.h"
#include "RotationMatrix.h"

registerMooseObject("MooseApp", InterfaceMaterialAvgBulkProp);

template <> InputParameters validParams<InterfaceMaterialAvgBulkProp>() {
  InputParameters params = validParams<Material>();
  params.set<bool>("is_interface_material") = true;
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "bulk_avg_mp_names", "names of average material proeprties");
  params.addRequiredParam<std::vector<UserObjectName>>("bulk_avg_mp_uo_names",
                                                       "names of averaging uo");
  params.addParam<std::vector<unsigned int>>(
      "stateful_flags",
      "a binary vector of 0 and 1  to declare if a material property need to "
      "be stafeul or not. 0 means non stateful (default), 1 means stateful");
  params.addClassDescription(
      "this material class is used in conjunction to a series of user objects "
      "to compute average material interface properties and add them as MP");

  return params;
}

InterfaceMaterialAvgBulkProp::InterfaceMaterialAvgBulkProp(
    const InputParameters &parameters)
    : Material(parameters),
      _bulk_avg_mp_names(
          getParam<std::vector<MaterialPropertyName>>("bulk_avg_mp_names")),
      _n_bulk_avg_mp(_bulk_avg_mp_names.size()),
      _bulk_avg_mp_uo_names(
          getParam<std::vector<UserObjectName>>("bulk_avg_mp_uo_names")),
      _n_bulk_avg_mp_uo_names(_bulk_avg_mp_uo_names.size()),
      _stateful_flags(
          parameters.isParamSetByUser("stateful_flags")
              ? getParam<std::vector<unsigned int>>("stateful_flags")
              : std::vector<unsigned int>(_n_bulk_avg_mp, 0)),
      _n_stateful_flags(_stateful_flags.size())

{

  if (_n_bulk_avg_mp != _n_bulk_avg_mp_uo_names)
    mooseError("InterfaceMaterialAvgBulkProp:: the number fo required"
               "averaged material properties does not "
               "match the number of provided averaging user objesct");

  if (_n_bulk_avg_mp != _n_stateful_flags)
    mooseError("InterfaceMaterialAvgBulkProp:: the number fo required"
               "averaged material properties, ",
               _n_bulk_avg_mp,
               ", does not "
               "match the number of provided stateful flags, ",
               _n_stateful_flags);
  // assgin user objects and declare associate material properties

  // assign user obejcts
  _bulk_avg_mp_uo.resize(_n_bulk_avg_mp);
  _avg_interface_mp.resize(_n_bulk_avg_mp);
  _avg_interface_mp_old.resize(_n_bulk_avg_mp);
  for (unsigned int i = 0; i < _n_bulk_avg_mp; ++i) {
    _bulk_avg_mp_uo[i] = &getUserObjectByName<ScalarBulkMPAcrossInterface_QP>(
        _bulk_avg_mp_uo_names[i]);
    _avg_interface_mp[i] = &declareProperty<Real>(_bulk_avg_mp_names[i]);
    if (_stateful_flags[i] == 1)
      _avg_interface_mp_old[i] =
          &getMaterialPropertyOldByName<Real>(_bulk_avg_mp_names[i]);
  }
}

void InterfaceMaterialAvgBulkProp::computeQpProperties() {
  // update average interface properties and save varaibles id

  for (unsigned int avg_mp_index = 0; avg_mp_index < _n_bulk_avg_mp;
       avg_mp_index++)
    (*_avg_interface_mp[avg_mp_index])[_qp] =
        _bulk_avg_mp_uo[avg_mp_index]->getAverage(_current_elem->id(),
                                                  _current_side, _qp);
}

void InterfaceMaterialAvgBulkProp::initQpStatefulProperties() {
  for (unsigned int avg_mp_index = 0; avg_mp_index < _n_bulk_avg_mp;
       avg_mp_index++)
    (*_avg_interface_mp[avg_mp_index])[_qp] = 0;
}
