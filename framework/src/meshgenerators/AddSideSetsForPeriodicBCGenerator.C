//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AddSideSetsForPeriodicBCGenerator.h"
#include "InputParameters.h"
#include "MooseTypes.h"
#include "MooseMeshUtils.h"
#include "CastUniquePointer.h"

#include "libmesh/remote_elem.h"

registerMooseObject("MooseApp", AddSideSetsForPeriodicBCGenerator);

defineLegacyParams(AddSideSetsForPeriodicBCGenerator);

InputParameters
AddSideSetsForPeriodicBCGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();

  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  params.addRequiredParam<RealVectorValue>(
      "dx_dy_dz", "A vector containg the shft of the periodic boudnaries in each direction.");
  params.addParam<Real>(
      "tolerance", 1e-10, "a fuzzy tolerance using when searching for periodic nodes");
  params.addClassDescription("MeshGenerator that creates sidesets for periodic BC");

  return params;
}

AddSideSetsForPeriodicBCGenerator::AddSideSetsForPeriodicBCGenerator(
    const InputParameters & parameters)
  : MeshGenerator(parameters),
    _input(getMesh("input")),
    _dx_dy_dz(getParam<RealVectorValue>("dx_dy_dz")),
    _tolerance(getParam<Real>("tolerance"))
{
}

bool
AddSideSetsForPeriodicBCGenerator::checkPointPeriodicity(const RealVectorValue & p1,
                                                         const RealVectorValue & p2,
                                                         const RealVectorValue & combination) const
{
  bool t = true;
  // RealVectorValue delta = (_dx_dy_dz * combination);
  for (unsigned int i = 0; i < 3; i++)
  {
    Real delta = _dx_dy_dz(i) * combination(i);
    if (std::abs(p1(i) + delta - p2(i)) > _tolerance)
    {
      t = false;
      break;
    }
  }
  return t;
}

std::vector<dof_id_type>
AddSideSetsForPeriodicBCGenerator::nodesOnSideGlobalId(const MeshBase & mesh,
                                                       const dof_id_type elem_id,
                                                       const unsigned int side) const
{
  const Elem * elem = mesh.elem_ptr(elem_id);
  const std::vector<unsigned int> local_nodes_on_side = elem->nodes_on_side(side);
  std::vector<unsigned> nodes_on_side;
  for (unsigned int i = 0; i < local_nodes_on_side.size(); i++)
    nodes_on_side.push_back(elem->node_id(local_nodes_on_side[i]));
  return nodes_on_side;
}

std::unique_ptr<MeshBase>
AddSideSetsForPeriodicBCGenerator::generate()
{
  std::unique_ptr<MeshBase> mesh = std::move(_input);
  const unsigned int dim = mesh->mesh_dimension();
  std::map<std::pair<dof_id_type, unsigned int>, std::vector<dof_id_type>> elem_face_node_map;
  std::unordered_set<dof_id_type> surface_nodes;
  for (const auto & elem : mesh->active_element_ptr_range())
  {
    for (unsigned int side = 0; side < elem->n_sides(); side++)
      if (elem->neighbor_ptr(side) == nullptr)
      {
        std::pair<dof_id_type, unsigned int> p = std::make_pair(elem->id(), side);
        const std::vector<dof_id_type> nodes_on_side_global =
            nodesOnSideGlobalId(*mesh, elem->id(), side);

        for (const auto & n_id : nodes_on_side_global)
          surface_nodes.insert(n_id);

        elem_face_node_map.insert(std::make_pair(p, nodes_on_side_global));
      }
  };
  // now we have all the free sides, we have all the nodes on the surface
  std::vector<RealVectorValue> combinations;
  combinations.push_back(Point(1, 0, 0));   // x
  combinations.push_back(Point(0, 1, 0));   // y
  combinations.push_back(Point(0, 0, 1));   // z
  combinations.push_back(Point(1, 1, 0));   // xy
  combinations.push_back(Point(1, -1, 0));  //-xy
  combinations.push_back(Point(1, 0, 1));   // xz
  combinations.push_back(Point(1, 0, -1));  //-xz
  combinations.push_back(Point(0, 1, 1));   // yz
  combinations.push_back(Point(0, 1, -1));  //-yz
  combinations.push_back(Point(1, 1, 1));   // xyz
  combinations.push_back(Point(1, -1, 1));  // x-yz
  combinations.push_back(Point(-1, 1, 1));  // -xyz
  combinations.push_back(Point(-1, -1, 1)); // -x-yz

  std::vector<std::pair<BoundaryName, BoundaryName>> sideset_names;
  sideset_names.push_back(std::make_pair<BoundaryName, BoundaryName>("per0_100", "per1_100"));
  sideset_names.push_back(std::make_pair<BoundaryName, BoundaryName>("per0_010", "per1_010"));
  sideset_names.push_back(std::make_pair<BoundaryName, BoundaryName>("per0_001", "per1_001"));
  sideset_names.push_back(std::make_pair<BoundaryName, BoundaryName>("per0_110", "per1_110"));
  sideset_names.push_back(std::make_pair<BoundaryName, BoundaryName>("per0_1-10", "per1_1-10"));
  sideset_names.push_back(std::make_pair<BoundaryName, BoundaryName>("per0_101", "per1_101"));
  sideset_names.push_back(std::make_pair<BoundaryName, BoundaryName>("per0_10-1", "per1_10-1"));
  sideset_names.push_back(std::make_pair<BoundaryName, BoundaryName>("per0_011", "per1_011"));
  sideset_names.push_back(std::make_pair<BoundaryName, BoundaryName>("per0_01-1", "per1_01-1"));
  sideset_names.push_back(std::make_pair<BoundaryName, BoundaryName>("per0_111", "per1_111"));
  sideset_names.push_back(std::make_pair<BoundaryName, BoundaryName>("per0_1-11", "per1_1-11"));
  sideset_names.push_back(std::make_pair<BoundaryName, BoundaryName>("per0_-111", "per1_-111"));
  sideset_names.push_back(std::make_pair<BoundaryName, BoundaryName>("per0_-1-11", "per1_-1-11"));
  // now we need to find which nodes pairs are periodic in the specified direction
  std::map<RealVectorValue, std::pair<BoundaryName, BoundaryName>> combination_sideset_pairs_map;
  for (unsigned int i = 0; i < combinations.size(); i++)
    combination_sideset_pairs_map[combinations[i]] =
        std::make_pair(sideset_names[i].first, sideset_names[i].second);

  std::map<BoundaryName, std::unordered_set<dof_id_type>> new_boundaries_map;
  // init the boudnary map
  for (const auto v : sideset_names)
  {
    new_boundaries_map[v.first] = std::unordered_set<dof_id_type>();
    new_boundaries_map[v.second] = std::unordered_set<dof_id_type>();
  }

  for (const auto & node_id_master : surface_nodes)
  {
    const Node * node_master = mesh->node_ptr(node_id_master);
    RealVectorValue xyz_master;
    xyz_master = (*node_master);

    for (const auto & node_id_slave : surface_nodes)
    {
      const Node * node_slave = mesh->node_ptr(node_id_slave);
      RealVectorValue xyz_slave;
      xyz_slave = (*node_slave);

      for (const auto & [c, ssn] : combination_sideset_pairs_map)
      {
        bool periodic = checkPointPeriodicity(xyz_master, xyz_slave, c);
        if (periodic)
        {
          auto it1 = new_boundaries_map.find(ssn.first);
          if (it1 != new_boundaries_map.end())
            it1->second.insert(node_id_master);

          auto it2 = new_boundaries_map.find(ssn.second);
          if (it2 != new_boundaries_map.end())
            it2->second.insert(node_id_slave);
        }
      }
    }
  }

  // now we have separated all the nodes in each direction, ideally we have nodesets
  std::vector<std::set<std::pair<dof_id_type, unsigned int>>> elem_side_list(dim * 2);

  std::map<BoundaryName, std::set<std::pair<dof_id_type, unsigned int>>> bname_elem_side_map;
  for (const auto v : sideset_names)
  {
    bname_elem_side_map[v.first] = std::set<std::pair<dof_id_type, unsigned int>>();
    bname_elem_side_map[v.second] = std::set<std::pair<dof_id_type, unsigned int>>();
  }

  // next step loop over the map and find which faces belong to which sideset
  for (const auto & [elem_face, nodes] : elem_face_node_map)
  {
    bool cant_find_face = true;
    for (const auto & [bname, nodes_in_set] : new_boundaries_map)
    {
      bool is_in = true;
      for (const auto & n : nodes)
      {
        is_in &= nodes_in_set.find(n) != nodes_in_set.end();
        if (!is_in)
          break;
      }

      if (is_in)
      {
        auto it1 = bname_elem_side_map.find(bname);
        if (it1 != bname_elem_side_map.end())
          it1->second.insert(elem_face);
        cant_find_face = false;
        break;
      }
    }
    if (cant_find_face)
      mooseError("AddSideSetsForPeriodicBCGenerator: can't find a periodic boundary for element " +
                 std::to_string(elem_face.first) + " face " + std::to_string(elem_face.second) +
                 " . Try to increase the tolerance for the search.");
  }

  // we are done seprating stuff along each direction, now we just need to create new sidesets.
  BoundaryInfo & boundary_info = mesh->get_boundary_info();
  const std::set<boundary_id_type> & currentBoundaryIds = boundary_info.get_boundary_ids();

  boundary_id_type max_bid = 0;
  for (const auto & b : currentBoundaryIds)
    if (b > max_bid)
      max_bid = b;
  max_bid += 1;

  for (const auto & [bname, elme_side_list] : bname_elem_side_map)
  {
    boundary_info.sideset_name(max_bid) = bname;
    for (auto & e_s : elme_side_list)
      boundary_info.add_side(e_s.first, e_s.second, max_bid);
    max_bid++;
  }

  return dynamic_pointer_cast<MeshBase>(mesh);
}
