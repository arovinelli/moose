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
class AssignFakeNeighborFromList;

template <>
InputParameters validParams<AssignFakeNeighborFromList>();

class AssignFakeNeighborFromList : public MeshGenerator
{
public:
  static InputParameters validParams();

  AssignFakeNeighborFromList(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  subdomain_id_type blockRestricteElementSubdomainID(const Elem * elem);
  void readFakeNeighborListFromFile();

  std::unique_ptr<MeshBase> & _input;
  const FileName _fake_neighbor_list_file_name;
  /// the total number of fake neighbors
  int _n_fake_neighbors = 0;

private:
  std::map<std::pair<dof_id_type, unsigned int>, std::pair<dof_id_type, unsigned int>>
      _fake_neighbor_map;
  std::map<std::pair<dof_id_type, unsigned int>, std::pair<dof_id_type, unsigned int>>
      _fake_neighbor_map_inverted;
};
