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

#ifndef BREAKMESHBYBLOCKBASE_H
#define BREAKMESHBYBLOCKBASE_H

#include "MeshModifier.h"

// forward declaration
class BreakMeshByBlockBase;

template <>
InputParameters validParams<BreakMeshByBlockBase>();

class BreakMeshByBlockBase : public MeshModifier
{
public:
  BreakMeshByBlockBase(const InputParameters & parameters);

  // method to override to implement other mesh splitting algorithms
  virtual void modify() override;

protected:
  /// the file_name from whence this mesh came
  std::string _file_name;
  /// the name of the new interface
  std::string _interface_name;
  /// the ID of the new interface
  BoundaryID _interface_id;
  /// the flag to split the interface by block
  bool _split_interface;

  /// check that if split_interface==true interface_id and interface_name are
  /// not set by the user. It also check that the provided interface_id is not
  /// already used
  void checkInputParameter();

  /// given the master and slave blocks this method return the appropriate
  /// boundary id and name
  void findBoundaryNameAndInd(const subdomain_id_type & /*masterBlockID*/,
                              const subdomain_id_type & /*slaveBlockID*/,
                              std::string & /*boundaryName*/,
                              BoundaryID & /*boundaryID*/,
                              BoundaryInfo & /*boundary_info*/);

  std::set<std::pair<std::string, BoundaryID>> _bName_bID_set;

private:
  /// this method finds the first free boundary id
  BoundaryID findFreeBoundaryId();

  /// this method generate the boundary name as
  /// "czm_bM_" + masterBlockID + "=>_bS_" + slaveBlockID;
  std::string generateBoundaryName(const subdomain_id_type & /*masterBlockID*/,
                                   const subdomain_id_type & /*slaveBlockID*/);

  /// this method save the boundary name/id pair
  void mapBoundaryIdAndBoundaryName(BoundaryID & /*boundaryID*/, std::string & /*boundaryName*/);
};

#endif // BREAKMESHBYBLOCKBASE_H
