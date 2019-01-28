//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MooseMesh.h"
#include "ScalarBulkMPAcrossInterface_QP.h"
registerMooseObject("MooseApp", ScalarBulkMPAcrossInterface_QP);

template <>
InputParameters
validParams<ScalarBulkMPAcrossInterface_QP>()
{
  InputParameters params = validParams<InterfaceUserObject>();
  params.addParam<MaterialPropertyName>("property_name",
                                        "the name of the scalar material property to averange "
                                        "across the interface");
  params.addParam<bool>("use_old_prop",
                        false,
                        "A Boolean to indicate whether the current or old "
                        "value of a material prop should be used.");
  return params;
}

ScalarBulkMPAcrossInterface_QP::ScalarBulkMPAcrossInterface_QP(const InputParameters & parameters)
  : InterfaceUserObject(parameters),
    _scalar_M(getParam<bool>("use_old_prop") ? getMaterialPropertyOld<Real>("property_name")
                                             : getMaterialProperty<Real>("property_name")),
    _scalar_S(getParam<bool>("use_old_prop") ? getNeighborMaterialPropertyOld<Real>("property_name")
                                             : getNeighborMaterialProperty<Real>("property_name"))

{
}

ScalarBulkMPAcrossInterface_QP::~ScalarBulkMPAcrossInterface_QP() {}

void
ScalarBulkMPAcrossInterface_QP::initialize()
{
  // define the boundary map nad retrieve element side and boundary_ID
  std::vector<std::tuple<dof_id_type, unsigned short int, boundary_id_type>> elem_side_bid =
      _mesh.buildSideList();

  // retrieve on which boudnary this UO operates
  std::set<BoundaryID> boundaryList = boundaryIDs();

  // clear map values
  _map_values.clear();

  // initialize the map_values looping over all the element and sides
  for (unsigned int i = 0; i < elem_side_bid.size(); i++)
  {
    // check if this boundary
    // if this element side is part of the boundary then add elements to the map
    if (boundaryList.find(std::get<2>(elem_side_bid[i])) != boundaryList.end())
    {

      // make pair
      std::pair<dof_id_type, unsigned int> elem_side_pair =
          std::make_pair(std::get<0>(elem_side_bid[i]), std::get<1>(elem_side_bid[i]));
      // initialize map elemenet
      std::vector<std::vector<Real>> var_values(0, std::vector<Real>(0));

      // add entry to the value map
      _map_values[elem_side_pair] = var_values;
    }
  }
}

void
ScalarBulkMPAcrossInterface_QP::execute()
{
  // find the entry on the map
  auto it = _map_values.find(std::make_pair(_current_elem->id(), _current_side));
  if (it != _map_values.end())
  {

    // insert two vector value for each qp
    auto & vec = _map_values[std::make_pair(_current_elem->id(), _current_side)];
    vec.resize(_qrule->n_points());

    // loop over qps and do stuff
    for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
    {
      vec[qp].resize(1);
      // compute average across interface
      vec[qp][0] = (_scalar_M[qp] + _scalar_S[qp]) / 2.;
    }
  }
  else
    mooseError("ScalarBulkMPAcrossInterface_QP::execute cannot fine the "
               "required element and side");
}

Real
ScalarBulkMPAcrossInterface_QP::getAverage(dof_id_type elem,
                                           unsigned int side,
                                           unsigned int qp) const
{
  auto value = _map_values.find(std::make_pair(elem, side));
  if (value != _map_values.end())
    return value->second[qp][0];
  else
    mooseError("ScalarBulkMPAcrossInterface_QP::getAverage can't find "
               "the given qp");
}
