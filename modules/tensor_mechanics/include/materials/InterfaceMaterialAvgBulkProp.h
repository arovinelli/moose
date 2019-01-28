//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef INTERFACEMATERIALAVGBULKPROP_H
#define INTERFACEMATERIALAVGBULKPROP_H

#include "Material.h"
#include "ScalarBulkMPAcrossInterface_QP.h"
class InterfaceMaterialAvgBulkProp;
template <> InputParameters validParams<InterfaceMaterialAvgBulkProp>();
/**
 *
 */
class InterfaceMaterialAvgBulkProp : public Material {
public:
  InterfaceMaterialAvgBulkProp(const InputParameters &parameters);

protected:
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;

  /// name of require averaged bulk mp
  std::vector<MaterialPropertyName> _bulk_avg_mp_names;
  /// number of require averaged bulk mp (inherited from _bulk_avg_mp_names
  /// length)
  const unsigned int _n_bulk_avg_mp;
  // /// userobjects returning the require bilk mp averaged
  std::vector<UserObjectName> _bulk_avg_mp_uo_names;
  const unsigned int _n_bulk_avg_mp_uo_names;

  //
  std::vector<unsigned int> _stateful_flags;
  const unsigned int _n_stateful_flags;
  std::vector<const ScalarBulkMPAcrossInterface_QP *> _bulk_avg_mp_uo;
  std::vector<MaterialProperty<Real> *> _avg_interface_mp;
  std::vector<const MaterialProperty<Real> *> _avg_interface_mp_old;
};

#endif // INTERFACEMATERIALAVGBULKPROP_H
