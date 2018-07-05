/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "BreakMeshByBlockManual_2DJunction.h"
#include "MooseMesh.h"

registerMooseObject("MooseApp", BreakMeshByBlockManual_2DJunction);

template <>
InputParameters
validParams<BreakMeshByBlockManual_2DJunction>()
{
  InputParameters params = validParams<BreakMeshByBlockManualBase>();
  params.addClassDescription("Manually split the mesh in 4ElementJunction.e."
                             "only works with REPLCIATED mesh");
  return params;
}

BreakMeshByBlockManual_2DJunction::BreakMeshByBlockManual_2DJunction(
    const InputParameters & parameters)
  : BreakMeshByBlockManualBase(parameters)
{
}

void
BreakMeshByBlockManual_2DJunction::modify()
{
  // if distributeb mesh raise an error
  _mesh_ptr->errorIfDistributedMesh(
      "BreakMeshByBlockManual_2DJunction only works on a REPLICATED mesh");

  checkInputParameter();

  updateElements();

  addInterfaceBoundary();
}

void
BreakMeshByBlockManual_2DJunction::updateElements()
{
  // specyfing nodes to duplciate

  dof_id_type curr_elem, local_node;

  curr_elem = 1;
  local_node = 0;
  duplicateAndSetLocalNode(curr_elem, local_node);

  curr_elem = 2;
  local_node = 0;
  duplicateAndSetLocalNode(curr_elem, local_node);

  curr_elem = 3;
  local_node = 0;
  duplicateAndSetLocalNode(curr_elem, local_node);

  curr_elem = 1;
  local_node = 1;
  duplicateAndSetLocalNode(curr_elem, local_node);

  curr_elem = 2;
  local_node = 2;
  duplicateAndSetLocalNode(curr_elem, local_node);

  curr_elem = 3;
  local_node = 1;
  duplicateAndSetLocalNode(curr_elem, local_node);

  curr_elem = 3;
  local_node = 2;
  duplicateAndSetLocalNode(curr_elem, local_node);
}

void
BreakMeshByBlockManual_2DJunction::addInterfaceBoundary()
{
  // construct boundary 100, which will contain the cohesive interface
  Elem * elem;
  BoundaryInfo & boundary_info = _mesh_ptr->getMesh().get_boundary_info();

  if (!_split_interface)
  {
    elem = _mesh_ptr->getMesh().elem(0);
    boundary_info.add_side(elem->id(), 0, _interface_id);
    boundary_info.add_side(elem->id(), 2, _interface_id);

    elem = _mesh_ptr->getMesh().elem(1);
    boundary_info.add_side(elem->id(), 2, _interface_id);

    elem = _mesh_ptr->getMesh().elem(2);
    boundary_info.add_side(elem->id(), 0, _interface_id);

    // rename the boundary
    boundary_info.sideset_name(_interface_id) = _interface_name;
  }
  else
  {
    elem = _mesh_ptr->getMesh().elem(0);
    boundary_info.add_side(elem->id(), 0, 0);
    boundary_info.sideset_name(0) = "bM1_bS2";

    boundary_info.add_side(elem->id(), 2, 5);
    boundary_info.sideset_name(5) = "bM1_bS3";

    elem = _mesh_ptr->getMesh().elem(1);
    boundary_info.add_side(elem->id(), 2, 6);
    boundary_info.sideset_name(6) = "bM2_bS4";

    elem = _mesh_ptr->getMesh().elem(2);
    boundary_info.add_side(elem->id(), 0, 7);
    boundary_info.sideset_name(7) = "bM3_bS4";
  }
}
