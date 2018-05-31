//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DisplacementJumpCohesiveInterface.h"
#include "MooseError.h"

registerMooseObject("TensorMechanicsApp", DisplacementJumpCohesiveInterface);

template <>
InputParameters
validParams<DisplacementJumpCohesiveInterface>()
{
  InputParameters params = validParams<InternalSideUserObject>();
  params += validParams<BoundaryRestrictableRequired>();
  // params.set<ExecFlagEnum>("execute_on") = EXEC_CUSTOM;
  // params.suppressParameter<ExecFlagEnum>("execute_on");
  params.addClassDescription("Compute the dispalcment jump across a cohesive interface");
  params.addRequiredCoupledVar("disp_x",
                               "variable containing the X component of"
                               "the dispalcement on the master side");
  params.addRequiredCoupledVar("disp_x_neighbor",
                               "variable containing the X component"
                               " of the dispalcement on the slave side");
  params.addCoupledVar("disp_y",
                       "variable containing the Y component of"
                       "the dispalcement on the master side");
  params.addCoupledVar("disp_y_neighbor",
                       "variable containing the Y "
                       "component of the dispalcement on the slave side");
  params.addCoupledVar("disp_z",
                       "variable containing the Z component of"
                       "the dispalcement on the master side");
  params.addCoupledVar("disp_z_neighbor",
                       "variable containing the Z "
                       "component of the dispalcement on the slave side");

  return params;
}

DisplacementJumpCohesiveInterface::DisplacementJumpCohesiveInterface(const InputParameters & params)
  : InternalSideUserObject(params),
    BoundaryRestrictableRequired(this, false), // false for applying to sidesets
    _disp_x(coupledValue("disp_x")),
    _disp_x_neighbor(coupledNeighborValue("disp_x_neighbor")),
    _disp_y(_mesh.dimension() >= 2 ? coupledValue("disp_y") : _zero),
    _disp_y_neighbor(_mesh.dimension() >= 2 ? coupledNeighborValue("disp_y_neighbor") : _zero),
    _disp_z(_mesh.dimension() >= 3 ? coupledValue("disp_z") : _zero),
    _disp_z_neighbor(_mesh.dimension() >= 3 ? coupledNeighborValue("disp_z_neighbor") : _zero)
{
}

RealVectorValue
DisplacementJumpCohesiveInterface::computeDisplacementJump()
{
  RealVectorValue Jump;
  Jump(0) = _disp_x_neighbor[_qp] - _disp_x[_qp];
  Jump(1) = _disp_y_neighbor[_qp] - _disp_y[_qp];
  Jump(2) = _disp_z_neighbor[_qp] - _disp_z[_qp];

  return Jump;
}

RealVectorValue
DisplacementJumpCohesiveInterface::RetrieveDisplacementJump(dof_id_type & elem, unsigned int & side)
{
  auto & vec = _map_values[std::pair<elem, side>];
  return vec;
}

// ovverride standard UO functions
void
DisplacementJumpCohesiveInterface::initialize()
{
  // find the number of qpoints for the element
  std::vector<dof_id_type> el;
  std::vector<unsigned short int> sl;
  std::vector<boundary_id_type> il;

  std::set<BoundaryID> boundaryList = boundaryIDs();
  _mesh.buildSideList(el, sl, il);
  _map_values.clear();

  for (unsigned int i = 0; i < il.size(); i++)
  {
    if (boundaryList.find(il[i]) != boundaryList.end())
    {
      std::pair<dof_id_type, unsigned int> elem_side_pair = std::make_pair(el[i], sl[i]);
      std::vector<RealVectorValue> var_values(0, RealVectorValue());
      _map_values.insert(std::make_pair(elem_side_pair, var_values));
    }
  }
}

void
DisplacementJumpCohesiveInterface::execute()
{

  auto & vec = _map_values[std::pair<_current_elem->id(), _current_side>];
  vec.resize(_qrule->npoints());
  vec[_qp] = computeDisplacementJump();

}

void
DisplacementJumpCohesiveInterface::finalize()
{
}

void
DisplacementJumpCohesiveInterface::threadJoin(const UserObject &)
{
}
