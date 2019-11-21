//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "JoinMeshByFaces.h"
#include "CastUniquePointer.h"

#include "libmesh/distributed_mesh.h"
#include "libmesh/elem.h"

#include <typeinfo>

registerMooseObject("MooseApp", JoinMeshByFaces);

template <>
InputParameters
validParams<JoinMeshByFaces>()
{
  InputParameters params = validParams<BreakMeshByBlockGeneratorBase>();
  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");
  params.addClassDescription("Break boundaries based on the subdomains to which their sides are "
                             "attached. Naming convention for the new boundaries will be the old "
                             "boundary name plus \"_to_\" plus the subdomain name. At the moment"
                             "this only works on REPLICATED mesh");
  return params;
}

JoinMeshByFaces::JoinMeshByFaces(const InputParameters & parameters)
  : BreakMeshByBlockGeneratorBase(parameters), _input(getMesh("input"))
{
  if (typeid(_input).name() == typeid(DistributedMesh).name())
    mooseError("JoinMeshByFaces only works with ReplicatedMesh.");
}

bool
JoinMeshByFaces::findOppositeElementFace(const Elem * elem,
                                         const std::vector<dof_id_type> opposite_nodes,
                                         unsigned int & opposite_elem_face) const
{
  bool found_opposite_face = false;
  for (unsigned int face_idx = 0; face_idx < elem->n_faces(); face_idx++)
    if (elem->neighbor_ptr(face_idx) == nullptr)
    {
      std::vector<dof_id_type> face_node_list(0);
      std::vector<dof_id_type> all_opposite_nodes(0);
      bool coincident_face = conincidentFace(elem, face_idx, face_node_list, all_opposite_nodes);
      std::cout << "Face face_idx:  " << face_idx << std::endl;
      if (coincident_face)
      {
        std::cout << "conincidentFace coincident_face:  " << face_idx << std::endl;
        ;
        std::cout << "findOppositeElementFace all_opposite_nodes ";
        for (auto n : all_opposite_nodes)
          std::cout << n << " ";
        std::cout << std::endl;

        std::cout << "findOppositeElementFace opposite_nodes ";
        for (auto n : opposite_nodes)
          std::cout << n << " ";
        std::cout << std::endl;

        std::vector<dof_id_type> node_intersection_list(opposite_nodes.size() +
                                                        all_opposite_nodes.size());
        auto it = std::set_intersection(opposite_nodes.begin(),
                                        opposite_nodes.end(),
                                        all_opposite_nodes.begin(),
                                        all_opposite_nodes.end(),
                                        node_intersection_list.begin());

        node_intersection_list.resize(it - node_intersection_list.begin());
        std::cout << "findOppositeElementFace node_intersection_list ";
        for (auto n : node_intersection_list)
          std::cout << n << " ";
        std::cout << std::endl;
        if (node_intersection_list.size() == opposite_nodes.size())
        {
          opposite_elem_face = face_idx;
          std::cout << "findOppositeElementFace opposite_elem_face " << opposite_elem_face
                    << "!!!!!!!!!!!!!" << std::endl;
          found_opposite_face = true;
          break;
        }
      }
    }
  return found_opposite_face;
}

std::unique_ptr<MeshBase>
JoinMeshByFaces::generate()
{
  std::unique_ptr<MeshBase> mesh = std::move(_input);

  // initialize the node to element map
  _node_to_elem_map.clear();
  for (const auto & elem : mesh->active_element_ptr_range())
    for (unsigned int n = 0; n < elem->n_nodes(); n++)
      _node_to_elem_map[elem->node_id(n)].push_back(elem->id());

  findCoincidentNodes(*mesh);
  findCoincidentFaces(*mesh);

  for (auto it : _coincident_faces_list)
  {
    Elem * elem = mesh->elem_ptr(it.first.first);
    elem->set_neighbor(it.first.second, mesh->elem_ptr(it.second.first));
  }

  for (auto it : _coincident_faces_list)
  {
    Elem * elem = mesh->elem_ptr(it.first.first);
    Elem * neigh = mesh->elem_ptr(it.second.first);
    std::cout << "elem " << elem->id() << "is the " << elem->which_neighbor_am_i(neigh)
              << "-th neighbor of element " << neigh->id() " face "
              << neigh->which_neighbor_am_i(elem) << std::endl;
    std::cout << "mapped data: elem " << it.first.first << " face " << it.first.second
              << " neighbor " << it.second.first << " face " << it.second.second << std::endl;
  }

  return dynamic_pointer_cast<MeshBase>(mesh);
}

void
JoinMeshByFaces::findCoincidentNodes(MeshBase & mesh)
{
  // clear map
  _coincident_nodes_map.clear();

  // loop over every ndoe pair to check equal coordinates
  for (const auto & node1 : mesh.local_node_ptr_range())
    for (const auto & node2 : mesh.local_node_ptr_range())
    {
      Point q;
      for (unsigned int k = 0; k < 3; ++k)
        q(k) = (*node2)(k);

      // check corridnates are ok and the node is not itself and add map element
      if (node1->absolute_fuzzy_equals(q) && (node1->id() != node2->id()))
      {
        _coincident_nodes_map[node1->id()].push_back(node2->id());
        std::cout << "node " << node1->id() << " == node " << node2->id() << std::endl;
      }
    }
}

bool
JoinMeshByFaces::findOppositeElement(const Elem * elem,
                                     const std::vector<dof_id_type> & face_nodes,
                                     dof_id_type & opposite_elem_id) const
{
  bool found_opposite_element = false;
  std::vector<dof_id_type> shared_elements(0);
  for (unsigned int i = 0; i < face_nodes.size(); i++)
  {
    std::vector<dof_id_type> opposite_nodes = _coincident_nodes_map.find(face_nodes[i])->second;
    std::vector<dof_id_type> opposite_elements(0);
    for (unsigned int j = 0; j < opposite_nodes.size(); j++)
    {
      auto it = _node_to_elem_map.find(opposite_nodes[j]);
      opposite_elements.insert(opposite_elements.end(), it->second.begin(), it->second.end());
    }
    std::sort(opposite_elements.begin(), opposite_elements.end());
    auto ip = std::unique(opposite_elements.begin(),
                          opposite_elements.begin() + opposite_elements.size());
    opposite_elements.resize(std::distance(opposite_elements.begin(), ip));

    if (shared_elements.size() == 0)
      shared_elements.assign(opposite_elements.begin(), opposite_elements.end());
    else
    {
      std::vector<dof_id_type> shared_elements_temp(shared_elements.size() +
                                                    opposite_elements.size());
      auto it = std::set_intersection(shared_elements.begin(),
                                      shared_elements.end(),
                                      opposite_elements.begin(),
                                      opposite_elements.end(),
                                      shared_elements_temp.begin());
      shared_elements_temp.resize(it - shared_elements_temp.begin());
      shared_elements.assign(shared_elements_temp.begin(), shared_elements_temp.end());
    }
  }
  if (shared_elements.size() == 1)
  {
    found_opposite_element = true;
    opposite_elem_id = shared_elements[0];
  }
  else
    mooseError("findOppositeElement: Can't fined the opposite lement, something is wrong");

  return found_opposite_element;
}

void
JoinMeshByFaces::findCoincidentFaces(MeshBase & mesh)
{
  // clear map
  _coincident_faces_list.clear();
  unsigned int n_coincident_face = 0;
  // loop over all the element and dine all the ones having a face with concident nodes using the
  // _node_to_elem_map map and the _coincident_nodes_map
  for (const auto & elem : mesh.active_element_ptr_range())
    for (unsigned int face_idx = 0; face_idx < elem->n_faces(); face_idx++)
    {
      std::vector<dof_id_type> face_node_list(0);
      std::vector<dof_id_type> all_duplicated_nodes(0);
      bool coincident_face = conincidentFace(elem, face_idx, face_node_list, all_duplicated_nodes);

      if (coincident_face)
      {
        std::cout << "elem " << elem->id() << " face " << face_idx << " is coincident. Nodes are ";
        for (auto n : face_node_list)
          std::cout << n << " ";
        std::cout << std::endl;
        n_coincident_face++;

        bool opposite_element_found = false;
        dof_id_type opposite_elem_id;
        opposite_element_found = findOppositeElement(elem, face_node_list, opposite_elem_id);

        if (opposite_element_found)
        {
          std::cout << " oppposite element id " << opposite_elem_id << std::endl;
          bool opposite_element_face_found = false;
          unsigned int opposite_element_face_id = 0;
          opposite_element_face_found = findOppositeElementFace(
              mesh.elem_ptr(opposite_elem_id), face_node_list, opposite_element_face_id);
          std::cout << "findCoincidentFaces opposite_element_face_found?"
                    << opposite_element_face_found << " ,face_id " << opposite_element_face_id
                    << "********************** " << std::endl;

          // save neighbor
          _coincident_faces_list[std::make_pair(elem->id(), face_idx)] =
              std::make_pair(opposite_elem_id, opposite_element_face_id);
        }
      }
    }
  std::cout << "I found " << n_coincident_face << "coincident faces" << std::endl;
}

bool
JoinMeshByFaces::conincidentFace(const Elem * elem,
                                 const unsigned int face_id,
                                 std::vector<dof_id_type> & face_node_list,
                                 std::vector<dof_id_type> & all_opposite_nodes) const
{
  bool coincident_face = true;
  face_node_list.resize(0);
  all_opposite_nodes.resize(0);
  for (auto local_n_id : elem->nodes_on_side(face_id))
  {
    auto it = _coincident_nodes_map.find(elem->node_id(local_n_id));
    if (it == _coincident_nodes_map.end())
    {
      coincident_face = false;
      face_node_list.resize(0);
      all_opposite_nodes.resize(0);
      break;
    }
    else
    {
      all_opposite_nodes.insert(all_opposite_nodes.end(), it->second.begin(), it->second.end());
      face_node_list.push_back(elem->node_id(local_n_id));
    }
  }
  if (coincident_face)
  {
    std::sort(face_node_list.begin(), face_node_list.end());
    std::sort(all_opposite_nodes.begin(), all_opposite_nodes.end());
    auto ip = std::unique(all_opposite_nodes.begin(), all_opposite_nodes.end());
    all_opposite_nodes.resize(std::distance(all_opposite_nodes.begin(), ip));
  }

  return coincident_face;
}

// void
// JoinMeshByFaces::addInterfaceBoundary(MeshBase & mesh)
// {
//   BoundaryInfo & boundary_info = mesh.get_boundary_info();
//
//   boundary_id_type boundaryID = findFreeBoundaryId(mesh);
//   std::string boundaryName = _interface_name;
//
//   // loop over boundary sides
//   for (auto & boundary_side_map : _new_boundary_sides_map)
//   {
//
//     // find the appropriate boundary name and id
//     //  given master and slave block ID
//     if (_split_interface)
//       findBoundaryNameAndInd(mesh,
//                              boundary_side_map.first.first,
//                              boundary_side_map.first.second,
//                              boundaryName,
//                              boundaryID,
//                              boundary_info);
//     else
//       boundary_info.sideset_name(boundaryID) = boundaryName;
//
//     // loop over all the side belonging to each pair and add it to the proper interface
//     for (auto & element_side : boundary_side_map.second)
//       boundary_info.add_side(element_side.first, element_side.second, boundaryID);
//   }
// }
