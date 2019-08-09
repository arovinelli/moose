//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterfaceQpValueUserObjectBase.h"
#include "MooseMesh.h"

template <>
InputParameters
validParams<InterfaceQpValueUserObjectBase>()
{
  InputParameters params = validParams<InterfaceValueUserObject>();
  params.addClassDescription("Test Interfae User Object computing and storing average values at "
                             "each QP across an interface");
  params.addParam<bool>("use_old_value", false, "each QP across an interface");
  return params;
}

InterfaceQpValueUserObjectBase::InterfaceQpValueUserObjectBase(const InputParameters & parameters)
  : InterfaceValueUserObject(parameters), _use_old_value(getParam<bool>("use_old_value"))
{
}

InterfaceQpValueUserObjectBase::~InterfaceQpValueUserObjectBase() {}

void
InterfaceQpValueUserObjectBase::initialize()
{

  _n_init += 1;

  if (_n_init == 1)
  {
    // define the boundary map and retrieve element side and boundary_ID
    std::vector<std::tuple<dof_id_type, unsigned short int, boundary_id_type>> elem_side_bid =
        _mesh.buildSideList();

    // retrieve on which boundary this UO operates
    std::set<BoundaryID> boundaryList = boundaryIDs();

    // clear map values
    _map_values.clear();

    // initialize the map_values looping over all the element and sides
    for (unsigned int i = 0; i < elem_side_bid.size(); i++)
    {
      // check if this element side is part of the boundary, if so add element side to the interface
      // map
      if (boundaryList.find(std::get<2>(elem_side_bid[i])) != boundaryList.end())
      {
        // make pair
        std::pair<dof_id_type, unsigned int> elem_side_pair =
            std::make_pair(std::get<0>(elem_side_bid[i]), std::get<1>(elem_side_bid[i]));
        // initialize map elemenet
        std::vector<std::vector<Real>> var_values(0, std::vector<Real>(0, 0));

        // add entry to the value map
        _map_values[elem_side_pair] = var_values;
      }
    }
  }
}

Real
InterfaceQpValueUserObjectBase::getQpValue(dof_id_type elem,
                                           unsigned int side,
                                           unsigned int qp) const
{
  auto data = _map_values.find(std::make_pair(elem, side));
  if (data != _map_values.end())
    return data->second[qp][0];
  else
    mooseError("InterfaceQpValueUserObjectBase::getQpValue can't find the given qp");
}
