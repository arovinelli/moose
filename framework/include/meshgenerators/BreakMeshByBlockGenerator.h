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
#include <unordered_set>

class BreakMeshByBlockGenerator;

template <>
InputParameters validParams<BreakMeshByBlockGenerator>();

class BreakMeshByBlockGenerator : public BreakMeshByBlockGeneratorBase
{
public:
  static InputParameters validParams();

  BreakMeshByBlockGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  subdomain_id_type blockRestricteElementSubdomainID(const Elem * elem);
  void writeFakeNeighborListToFile(MeshBase & mesh) const;

  std::unique_ptr<MeshBase> & _input;
  std::vector<SubdomainID> _block;
  std::unordered_set<SubdomainID> _block_set;
  const bool _block_restricted;
  const bool _add_transition_interface;
  const bool _split_transition_interface;
  const bool _write_fake_neighbor_list_to_file;
  const FileName _fake_neighbor_list_file_name;
  /// the total number of fake neighbors
  int _n_fake_neighbors = 0;

private:
  /// generate the new boundary interface
  void addInterfaceBoundary(MeshBase & mesh);

  std::set<std::pair<subdomain_id_type, subdomain_id_type>> _neighboring_block_list;
  std::map<std::pair<subdomain_id_type, subdomain_id_type>,
           std::set<std::pair<dof_id_type, unsigned int>>>
      _new_boundary_sides_map;
};
