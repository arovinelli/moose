//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CohesiveInterfaceMaterialBase.h"
#include "Assembly.h"

#include "RotationMatrix.h"

template <>
InputParameters
validParams<CohesiveInterfaceMaterialBase>()
{
  InputParameters params = validParams<Material>();

  params.addRequiredParam<UserObjectName>(
      "uo_neighbor_bulk_properties_and_variables",
      "the name of the user object retrieveing neighboring bulk material properties and variables");
  params.addParam<bool>("is_interface_material",
                        true,
                        "this parameter shouls always be TRUE and never changed by the User");
  params.addClassDescription(
      "this material class is the base class from which a cohesive zone"
      "material shoul always be derived. It contains kinematic functions"
      "to convert the disaplcement jump in from global to local coordinates");
  return params;
}

CohesiveInterfaceMaterialBase::CohesiveInterfaceMaterialBase(const InputParameters & parameters)
  : Material(parameters),
    _normals(_assembly.normals()),
    _uo_neighbor_bulk_properties_and_variables(
        getUserObject<DispJumpAcrossInterface>("uo_neighbor_bulk_properties_and_variables"))

{
}

void
CohesiveInterfaceMaterialBase::computeQpProperties()
{
  mooseError("CohesiveInterfaceMaterialBase::computeQpProperties should never be called."
             "Only its derived classes are allowed to do so");
}

void
CohesiveInterfaceMaterialBase::initQpStatefulProperties()
{

  mooseError("CohesiveInterfaceMaterialBase::initQpStatefulProperties should never be called."
             "Only its derived classes are allowed to do so");
}
