//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DispJumpAcrossInterface.h"
#include "MooseError.h"
// #include "Assembly.h"
// #include "RotationMatrix.h"

registerMooseObject("TensorMechanicsApp", DispJumpAcrossInterface);

template <>
InputParameters
validParams<DispJumpAcrossInterface>()
{
  InputParameters params = validParams<InterfaceUserObject>();

  params.addClassDescription("Compute the dispalcment jump across an interface");
  params.addRequiredCoupledVar("disp_x",
                               "variable containing the X component of"
                               "the displacement on the master side");
  params.addRequiredCoupledVar("disp_x_neighbor",
                               "variable containing the X component"
                               " of the displacement on the slave side");
  params.addCoupledVar("disp_y",
                       "variable containing the Y component of"
                       "the displacement on the master side");
  params.addCoupledVar("disp_y_neighbor",
                       "variable containing the Y "
                       "component of the displacement on the slave side");
  params.addCoupledVar("disp_z",
                       "variable containing the Z component of"
                       "the displacement on the master side");
  params.addCoupledVar("disp_z_neighbor",
                       "variable containing the Z "
                       "component of the displacement on the slave side");

  return params;
}

DispJumpAcrossInterface::DispJumpAcrossInterface(const InputParameters & params)
  : InterfaceUserObject(params),
    _disp_x(coupledValue("disp_x")),
    _disp_x_neighbor(coupledNeighborValue("disp_x_neighbor")),
    _disp_y(_mesh.dimension() >= 2 ? coupledValue("disp_y") : _zero),
    _disp_y_neighbor(_mesh.dimension() >= 2 ? coupledNeighborValue("disp_y_neighbor") : _zero),
    _disp_z(_mesh.dimension() >= 3 ? coupledValue("disp_z") : _zero),
    _disp_z_neighbor(_mesh.dimension() >= 3 ? coupledNeighborValue("disp_z_neighbor") : _zero)
{
}

DispJumpAcrossInterface::~DispJumpAcrossInterface() {}

// ovverride standard UO functions
void
DispJumpAcrossInterface::initialize()
{

  // get the boundary list on which this user object works
  std::set<BoundaryID> boundaryList = boundaryIDs();

  // _mesh.buildSideList(el, sl, il);
  // initialize and fill element side boundary id lists
  std::vector<std::tuple<dof_id_type, unsigned short int, boundary_id_type>> elem_side_bid;
  elem_side_bid = _mesh.buildSideList();

  // reset map values
  _map_values.clear();

  // initialize the map in which values to return are stored
  for (unsigned int i = 0; i < elem_side_bid.size(); i++)
  {
    if (boundaryList.find(std::get<2>(elem_side_bid[i])) != boundaryList.end())
    {
      std::pair<dof_id_type, unsigned int> elem_side_pair =
          std::make_pair(std::get<0>(elem_side_bid[i]), std::get<1>(elem_side_bid[i]));
      std::vector<std::vector<RealVectorValue>> var_values(0, std::vector<RealVectorValue>(0));
      _map_values[elem_side_pair] = var_values;
    }
  }
}

void
DispJumpAcrossInterface::execute()
{
  // this method is executed for each element therefore we need to loop over _qp
  auto it = _map_values.find(std::make_pair(_current_elem->id(), _current_side));
  if (it != _map_values.end())
  {

    // retrieve the propoer locala in the map
    auto & vec = _map_values[std::make_pair(_current_elem->id(), _current_side)];
    // resize the temporary vector accoridng to the number of qp
    vec.resize(_qrule->n_points());

    // cycle over qp
    for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
    {
      // resize container jsut to be sure
      vec[qp].resize(1, RealVectorValue());

      // store dispalcment jump
      vec[qp][0] = computeDisplacementJump(qp);
      // // store normal at qp on the master interface
      // vec[qp][1] = _normals[qp];
    }
  }
  else
    mooseError("DispJumpAcrossInterface::execute can't find the saved value");
}

void
DispJumpAcrossInterface::threadJoin(const UserObject &)
{
}

void
DispJumpAcrossInterface::finalize()
{
}

RealVectorValue
DispJumpAcrossInterface::computeDisplacementJump(unsigned int qp)
{

  RealVectorValue Jump;
  // compute the displacement jump
  Jump(0) = _disp_x_neighbor[qp] - _disp_x[qp];
  Jump(1) = _disp_y_neighbor[qp] - _disp_y[qp];
  Jump(2) = _disp_z_neighbor[qp] - _disp_z[qp];

  return Jump;
}

RealVectorValue
DispJumpAcrossInterface::getDisplacementJump(const dof_id_type elem_id,
                                             const dof_id_type side_id,
                                             const unsigned int qp)
{
  // find the appropriate map acces
  auto it = _map_values.find(std::make_pair(elem_id, side_id));
  if (it != _map_values.end())
    // return the displacement jump (sotred as 0 element)
    return it->second[qp][0];
  else
    mooseError("DispJumpAcrossInterface::getDisplacementJump can't fidn the saved value");
}
