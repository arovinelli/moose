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

class JoinMeshByFaces;

template <>
InputParameters validParams<JoinMeshByFaces>();

class JoinMeshByFaces : public BreakMeshByBlockGeneratorBase
{
public:
  JoinMeshByFaces(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  std::unique_ptr<MeshBase> & _input;
  std::unordered_map<dof_id_type, std::vector<dof_id_type>> _coincident_nodes_map;
  std::unordered_map<std::pair<dof_id_type, unsigned int>, std::pair<dof_id_type, unsigned int>>
      _coincident_faces_list;
  std::unordered_map<dof_id_type, std::vector<dof_id_type>> _node_to_elem_map;

  void findCoincidentNodes(MeshBase & mesh);
  void findCoincidentFaces(MeshBase & mesh);
  void findCoincidentElementFaces(dof_id_type elem, std::vector<dof_id_type> opposite_nodes);
  bool conincidentFace(const Elem * elem,
                       const unsigned int face_id,
                       std::vector<dof_id_type> & face_node_list,
                       std::vector<dof_id_type> & all_opposite_nodes) const;

  bool findOppositeElementFace(const Elem * elem,
                               const std::vector<dof_id_type> opposite_nodes,
                               unsigned int & opposite_elem_face) const;

  bool findOppositeElement(const Elem * elem,
                           const std::vector<dof_id_type> & face_nodes,
                           dof_id_type & opposite_elem_id) const;

private:
  // /// generate the new boundary interface
  // void addInterfaceBoundary(MeshBase & mesh);

  // std::set<std::pair<subdomain_id_type, subdomain_id_type>> _neighboring_block_list;
  // std::map<std::pair<subdomain_id_type, subdomain_id_type>,
  //          std::set<std::pair<dof_id_type, unsigned int>>>
  //     _new_boundary_sides_map;
};
