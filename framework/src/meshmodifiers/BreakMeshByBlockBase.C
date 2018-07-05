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

#include "BreakMeshByBlockBase.h"
#include "MooseMesh.h"

template <>
InputParameters
validParams<BreakMeshByBlockBase>()
{
  InputParameters params = validParams<MeshModifier>();
  params.addParam<std::string>("interface_name",
                               "interface",
                               "the name of"
                               "the new interface");
  params.addParam<BoundaryID>("interface_id",
                              100,
                              "the ID of the new"
                              "interface");
  params.addParam<bool>("split_interface",
                        false,
                        "If true, it create a "
                        "different interface for each block pair");
  params.addClassDescription("This is the base class used to split a monolithic"
                             "mesh by blocks pairs");
  return params;
}

BreakMeshByBlockBase::BreakMeshByBlockBase(const InputParameters & parameters)
  : MeshModifier(parameters),
    _interface_name(getParam<std::string>("interface_name")),
    _interface_id(getParam<BoundaryID>("interface_id")),
    _split_interface(getParam<bool>("split_interface"))
{
}

void
BreakMeshByBlockBase::modify()
{
  mooseError("BreakMeshByBlockBase should never be called directly!"
             "Always use one of its the derived classes");
}

void
BreakMeshByBlockBase::checkInputParameter()
{
  // check input consistency
  if (getParam<bool>("split_interface") &&
      (_pars.isParamSetByUser("interface_id") || _pars.isParamSetByUser("interface_name")))
  {
    mooseError("if split_interface == true, interface_name and/or interface_id"
               " cannot be specified by the user");
  }
  else
  {
    // check that the provided interface_id is not already used
    const std::set<BoundaryID> & currentBoundaryIds = _mesh_ptr->getBoundaryIDs();
    if (currentBoundaryIds.count(getParam<BoundaryID>("interface_id")) != 0)
      mooseError("BreakMeshByBlockBase::checkInputParameter. "
                 "A boundary with the same interface_id already exists "
                 "in the mesh. Please select a different interface_id.");
  }
}

BoundaryID
BreakMeshByBlockBase::findFreeBoundaryId()
{
  const std::set<BoundaryID> & currentBoundaryIds = _mesh_ptr->getBoundaryIDs();
  bool freeBoundaryNotFound = true;
  BoundaryID freeId;
  for (freeId = 0; freeId < std::numeric_limits<BoundaryID>::max(); freeId++)
  {
    if (currentBoundaryIds.count(freeId) == 0)
    {
      // bid is not in the set, boundaryID is free
      freeBoundaryNotFound = false;
      break;
    }
  }

  if (freeBoundaryNotFound)
    mooseError("Too many boundaries. Maximum limit exceeded!");

  return freeId;
}

std::string
BreakMeshByBlockBase::generateBoundaryName(const subdomain_id_type & masterBlockID,
                                           const subdomain_id_type & slaveBlockID)
{
  return "bM" + std::to_string(masterBlockID) + "_bS" + std::to_string(slaveBlockID);
}

void
BreakMeshByBlockBase::mapBoundaryIdAndBoundaryName(BoundaryID & boundaryID,
                                                   std::string & boundaryName)
{
  _bName_bID_set.insert(std::pair<std::string, int>(boundaryName, boundaryID));
}

void
BreakMeshByBlockBase::findBoundaryNameAndInd(const subdomain_id_type & masterBlockID,
                                             const subdomain_id_type & slaveBlockID,
                                             std::string & boundaryName,
                                             BoundaryID & boundaryID,
                                             BoundaryInfo & boundary_info)
{
  // TODO this method should be multiprocessor and multithread safe
  // rememeber to check with the developer team

  // mpi barrier
  // first check which boundary name will be created
  boundaryName = generateBoundaryName(masterBlockID, slaveBlockID);

  // check if the boundary name already exist
  bool checkBoundaryAlreadyExist = false;
  for (auto b : _bName_bID_set)
  {
    if (b.first.compare(boundaryName) == 0)
    {
      boundaryID = b.second;
      checkBoundaryAlreadyExist = true;
    }
  }

  if (checkBoundaryAlreadyExist)
  {
    // mpi barrier end
    return;
  }
  else
  {
    boundaryID = findFreeBoundaryId();
    mapBoundaryIdAndBoundaryName(boundaryID, boundaryName);

    // rename the boundary
    boundary_info.sideset_name(boundaryID) = boundaryName;

    // mpi broadcast boundary info and _bName_bID_set

    // mpi barrier end
    return;
  }
}
