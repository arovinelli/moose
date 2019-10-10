//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseObject.h"
#include "Restartable.h"

// Included so mesh generators don't need to include this when constructing MeshBase objects
#include "MooseMesh.h"

#include "libmesh/mesh_base.h"
#include "libmesh/parameters.h"

// Forward declarations
class MeshGenerator;
class MooseMesh;

template <>
InputParameters validParams<MeshGenerator>();

/**
 * MeshGenerators are objects that can modify or add to an existing mesh.
 */
class MeshGenerator : public MooseObject, public MeshMetaDataInterface, public Restartable
{
public:
  /**
   * Constructor
   *
   * @param parameters The parameters object holding data for the class to use.
   */
  MeshGenerator(const InputParameters & parameters);

  /**
   * Generate / modify the mesh
   *
   */
  virtual std::unique_ptr<MeshBase> generate() = 0;

  /**
   * Return the MeshGenerators that must run before this MeshGenerator
   */
  std::vector<std::string> & getDependencies() { return _depends_on; }

protected:
  template <typename T>
  void setProperty(const std::string & name, const T & value);

  /**
   * Takes the name of a MeshGeneratorName parameter and then gets a pointer to the
   * Mesh that MeshGenerator is going to create.
   *
   * NOTE: You MUST catch this by reference!
   *
   * @return The Mesh generated by that MeshGenerator
   */
  std::unique_ptr<MeshBase> & getMesh(const std::string & input_mesh_generator_parameter_name);

  /**
   * Takes the name of another MeshGenerator directly.
   *
   * NOTE: You MUST catch this by reference!
   *
   * @return The Mesh generated by that MeshGenerator
   */
  std::unique_ptr<MeshBase> &
  getMeshByName(const MeshGeneratorName & input_mesh_generator_parameter_name);

  /// References to the mesh and displaced mesh (currently in the ActionWarehouse)
  std::shared_ptr<MooseMesh> & _mesh;

private:
  /// A list of generators that are required to run before this generator may run
  std::vector<std::string> _depends_on;

  /// A nullptr to use for when inputs aren't specified
  std::unique_ptr<MeshBase> _null_mesh = nullptr;

  Parameters & _mesh_meta_data;
};

template <typename T>
void
MeshGenerator::setProperty(const std::string & name, const T & value)
{
  _mesh_meta_data.set<T>(name) = value;
}
