//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "BreakMeshByBlockGeneratorBase.h"

class ParallelSyncExample;

template <>
InputParameters validParams<ParallelSyncExample>();

class ParallelSyncExample : public MeshGenerator
{
public:
  static InputParameters validParams();

  ParallelSyncExample(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  std::unique_ptr<MeshBase> & _input;

private:
  std::vector<dof_id_type>
  findNodeMultiplicity(const std::map<dof_id_type, std::vector<dof_id_type>> & node_to_elem_map,
                       const MeshBase & mesh);

  void getMyDataParallel(const std::vector<dof_id_type> & node_multiplicity_local,
                         const MeshBase & mesh);
};
