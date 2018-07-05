//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BreakMeshByBlock.h"
#include "MooseMesh.h"

registerMooseObject("MooseApp", BreakMeshByBlock);

template <>
InputParameters
validParams<BreakMeshByBlock>()
{
  InputParameters params = validParams<BreakMeshByBlockBase>();
  params.addClassDescription("Break boundaries based on the subdomains to which their sides are "
                             "attached. Naming convention for the new boundaries will be the old "
                             "boundary name plus \"_to_\" plus the subdomain name. At the moment"
                             "this only works on REPLICATED mesh");
  params.addParam<std::vector<BoundaryName>>(
      "boundaries", "Boundaries to be broken. Default means to break all boundaries");
  return params;
}

BreakMeshByBlock::BreakMeshByBlock(const InputParameters & parameters)
  : BreakMeshByBlockBase(parameters)
{
}

void
BreakMeshByBlock::modify()
{
  // TODO remove when distributed MESH capabilities are implemented
  _mesh_ptr->errorIfDistributedMesh("BreakMeshByBlock only works on a REPLICATED mesh");

  checkInputParameter();
  buildNodeSupport();
  buildInterfacialNodes();
  duplicateNodes();
  tearElements();
  addInterfaceBoundary();
}

void
BreakMeshByBlock::buildNodeSupport()
{
  for (MeshBase::element_iterator it = _mesh_ptr->getMesh().elements_begin();
       it != _mesh_ptr->getMesh().elements_end();
       ++it)
  {
    Elem * elem = *it;
    for (unsigned int i = 0; i < elem->n_nodes(); ++i)
      _node_support[elem->node_id(i)].push_back(elem->id());
  }
}

void
BreakMeshByBlock::buildInterfacialNodes()
{

  std::set<subdomain_id_type> matSet;

  // loop over all the nodes
  for (MeshBase::node_iterator it = _mesh_ptr->getMesh().nodes_begin();
       it != _mesh_ptr->getMesh().nodes_end();
       ++it)
  {
    Node * node = *it;
    std::vector<dof_id_type> support = _node_support[node->id()];
    unsigned int suppCount = support.size();

    if (suppCount == 1)
      continue;

    for (unsigned int ie = 0; ie < suppCount; ie++)
    {
      subdomain_id_type imat = _mesh_ptr->getMesh().elem(support[ie])->subdomain_id();
      matSet.insert(imat);
    }

    unsigned int matCount = matSet.size();

    // if interface elements are created everywhere then
    // duplicity = nodal support. Else duplicity = materials.
    unsigned int duplicity = matCount;

    if (matCount != 1) // interfacial node
    {
      _node_duplicity[node->id()] = duplicity;

      _duplicated_node_materials[node->id()].insert(
          _duplicated_node_materials[node->id()].end(), matSet.begin(), matSet.end());
    }

    matSet.clear(); // clear for the next node
  }
}

void
BreakMeshByBlock::duplicateNodes()
{
  int idd(0);

  // after we have found the node duplicity we duplicate all the nodes as needed
  for (auto & dn : _node_duplicity)
  {
    _duplicated_node[dn.first].push_back(dn.first);

    // duplicate the nodes dn.second - 1 times (the original node remin assigned to one element)
    for (unsigned int id = 0; id < dn.second - 1; id++)
    {
      idd++;

      Node * new_node =
          Node::build(_mesh_ptr->getMesh().node(dn.first), _mesh_ptr->getMesh().n_nodes())
              .release();
      new_node->processor_id() = _mesh_ptr->getMesh().node(dn.first).processor_id();
      _mesh_ptr->getMesh().add_node(new_node);

      _duplicated_node[dn.first].push_back(new_node->id());
    }
  }
}

void
BreakMeshByBlock::tearElements()
{
  for (auto & dn : _duplicated_node) // loop over diplicated nodes
  {
    std::vector<dof_id_type> support = _node_support[dn.first];
    unsigned int suppCount = support.size();

    Elem * ref_elem = _mesh_ptr->getMesh().elem(support[0]);

    // block id of the first element in the suppor of the dn
    subdomain_id_type ref_mat = ref_elem->subdomain_id();

    for (unsigned int ie = 1; ie < suppCount; ie++) // loop over adjacent elments in the support
    {

      // 1 -----> get the currenet element
      Elem * i_elem = _mesh_ptr->getMesh().elem(support[ie]);
      subdomain_id_type imat = i_elem->subdomain_id();

      if (imat == ref_mat)
        continue;
      // if material is different then we ne need to duplicate the node and store the boundary
      // pair

      // 1 -----> store boundary pair by block ID
      std::pair<subdomain_id_type, subdomain_id_type> materials_pair;

      if (imat > ref_mat)
        materials_pair = std::make_pair(ref_mat, imat);
      else
        materials_pair = std::make_pair(imat, ref_mat);

      // store the boundary pair
      _boundary_pairs.insert(materials_pair);

      // 2 -----> find the shared node and assign it to the appropriate duplicate node
      unsigned int imat_id = 0;
      for (unsigned int i = 0; i < i_elem->n_nodes(); ++i)
      {

        if (i_elem->node_id(i) == dn.first) // dn.first is the index of the original duplicated node
        {

          std::vector<subdomain_id_type> mats = _duplicated_node_materials[dn.first];
          // find the appropriate duplicate node by cycling materials
          for (unsigned int i = 0; i < mats.size(); i++)
          {
            if (imat == mats[i])
            {
              imat_id = i;
              break;
            }
          }

          // assign new node to the element
          i_elem->set_node(i) = _mesh_ptr->getMesh().node_ptr((dn.second)[imat_id]);
        }
      }

      // 3 -----> find adjacent faces (if any) always sorted according to the material id
      // element on the master side of the boundary (i.e. the one that
      // will contain the interface)
      Elem * boundary_elem;

      // element on slave side of the boundary
      Elem * adjacent_elem;

      // the material of the boundary side
      // subdomain_id_type boundary_mat;

      // identify master(boundary) and slave(adjcent) element by looking at the
      // material ID the master side is always the one with the lowest ID
      if (imat > ref_mat)
      {
        boundary_elem = ref_elem;
        adjacent_elem = i_elem;
        // boundary_mat = ref_mat;
      }
      else
      {
        boundary_elem = i_elem;
        adjacent_elem = ref_elem;
        // boundary_mat = imat;
      }

      // save all boundary faces
      for (unsigned int i = 0; i < boundary_elem->n_sides(); i++)
      { // loop over faces
        if (boundary_elem->neighbor(i) != nullptr)
        {
          if (boundary_elem->neighbor(i)->id() == adjacent_elem->id())
          {
            _boundary_sides[materials_pair].insert(std::make_pair(boundary_elem->id(), i));
            break; // between two elements only one face can be common
          }
        }
      }
    }
  }
}

void
BreakMeshByBlock::addInterfaceBoundary()
{
  BoundaryInfo & boundary_info = _mesh_ptr->getMesh().get_boundary_info();

  BoundaryID boundaryID = _interface_id;
  std::string boundaryName = _interface_name;

  // loop over boundary sides
  for (auto & bs : _boundary_sides)
  {

    if (_split_interface)
    {
      // find the appropriate boundary name and id
      //  given master and slave block ID
      findBoundaryNameAndInd(
          bs.first.first, bs.first.second, boundaryName, boundaryID, boundary_info);
    }
    else
      boundary_info.sideset_name(boundaryID) = boundaryName;

    // loop over all the side belonging to each pair and add it to the proper interface
    for (auto & es : bs.second)
    {
      boundary_info.add_side(es.first, es.second, boundaryID);
    }
  }
}
