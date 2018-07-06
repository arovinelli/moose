//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef COHESIVEINTERFACEMATERIALBASE_H
#define COHESIVEINTERFACEMATERIALBASE_H

#include "Material.h"
#include "DispJumpAcrossInterface.h"

class CohesiveInterfaceMaterialBase;

template <>
InputParameters validParams<CohesiveInterfaceMaterialBase>();

/// base class for implementing cohesive zone materials
class CohesiveInterfaceMaterialBase : public Material
{
public:
  CohesiveInterfaceMaterialBase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;

  /// normal of the current side of the interface
  const MooseArray<Point> & _normals;

  /// the user object grabbing all the variables and stateful material properties
  /// required to update the state of this material
  const DispJumpAcrossInterface & _uo_neighbor_bulk_properties_and_variables;

  /// teh variable representing the displacement jump across the interface
  RealVectorValue _displacement_jump;
};

#endif // COHESIVEINTERFACEMATERIALBASE_H
