//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MeshGenerator.h"

// Forward declarations
class AddSideSetsForPeriodicBCGenerator;

template <>
InputParameters validParams<AddSideSetsForPeriodicBCGenerator>();

/**
 * MeshGenerator that creates a sideset composed of the nodes located between
 * two or more subdomains.
 */
class AddSideSetsForPeriodicBCGenerator : public MeshGenerator
{
public:
  static InputParameters validParams();

  AddSideSetsForPeriodicBCGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  std::unique_ptr<MeshBase> & _input;
  std::vector<dof_id_type> nodesOnSideGlobalId(const MeshBase & /*mesh*/,
                                               const dof_id_type /*elem_id*/,
                                               const unsigned int /*side*/) const;

  bool checkPointPeriodicity(const RealVectorValue & /*p1*/,
                             const RealVectorValue & /*p2*/,
                             const RealVectorValue & /*combination*/) const;

  const RealVectorValue _dx_dy_dz;
  const Real _tolerance;
};
