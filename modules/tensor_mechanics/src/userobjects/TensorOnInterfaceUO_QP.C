//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TensorOnInterfaceUO_QP.h"
#include "MooseMesh.h"
registerMooseObject("TensorMechanicsApp", TensorOnInterfaceUO_QP);

template <>
InputParameters
validParams<TensorOnInterfaceUO_QP>()
{
  InputParameters params = validParams<InterfaceUserObject>();
  params.addRequiredParam<MaterialPropertyName>(
      "tensor_mp_name", "the name of the tensor we want to averange along the interface");
  params.addParam<bool>(
      "use_old_prop",
      false,
      "A Boolean to indicate whether the current or old value of a material prop should be used.");
  return params;
}

TensorOnInterfaceUO_QP::TensorOnInterfaceUO_QP(const InputParameters & parameters)
  : InterfaceUserObject(parameters),
    _tensor_M(getParam<bool>("use_old_prop")
                  ? getMaterialPropertyOld<RankTwoTensor>("tensor_mp_name")
                  : getMaterialProperty<RankTwoTensor>("tensor_mp_name")),
    _tensor_S(getParam<bool>("use_old_prop")
                  ? getNeighborMaterialPropertyOld<RankTwoTensor>("tensor_mp_name")
                  : getNeighborMaterialProperty<RankTwoTensor>("tensor_mp_name"))

{
}

TensorOnInterfaceUO_QP::~TensorOnInterfaceUO_QP() {}

void
TensorOnInterfaceUO_QP::initialize()
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
      std::vector<std::vector<RankTwoTensor>> var_values(0, std::vector<RankTwoTensor>(0));

      // add entry to the value map
      _map_values[elem_side_pair] = var_values;
    }
  }
}

void
TensorOnInterfaceUO_QP::execute()
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
      vec[qp].resize(2);
      // compute tensor average
      vec[qp][0] = _tensor_M[qp];
      // compute tensor difference
      vec[qp][1] = _tensor_S[qp];
    }
  }
  else
    mooseError("TensorOnInterfaceUO_QP::execute cannot fine the required element and side");
}

RankTwoTensor
TensorOnInterfaceUO_QP::getTensorMaster(dof_id_type elem, unsigned int side, unsigned int qp) const
{
  auto Tensor = _map_values.find(std::make_pair(elem, side));
  if (Tensor != _map_values.end())
    return Tensor->second[qp][0];
  else
    mooseError("TensorOnInterfaceUO_QP::getTensorAverage can't find the given qp");
}

RankTwoTensor
TensorOnInterfaceUO_QP::getTensorSlave(dof_id_type elem, unsigned int side, unsigned int qp) const
{
  auto Tensor = _map_values.find(std::make_pair(elem, side));
  if (Tensor != _map_values.end())
    return Tensor->second[qp][1];
  else
    mooseError("TensorOnInterfaceUO_QP::getTensorDifference can't find the given qp");
}
