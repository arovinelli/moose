//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterfaceQpValueUserObject.h"
#include "MooseMesh.h"
registerMooseObject("MooseApp", InterfaceQpValueUserObject);

template <>
InputParameters
validParams<InterfaceQpValueUserObject>()
{
  InputParameters params = validParams<InterfaceQpValueUserObjectBase>();
  params.addRequiredCoupledVar("var", "The variable name");
  params.addCoupledVar("var_neighbor", "The variable name");
  params.addClassDescription("Test Interfae User Object computing and storing average values at "
                             "each QP across an interface");
  params.addParam<bool>("use_old_value", false, "each QP across an interface");
  return params;
}

InterfaceQpValueUserObject::InterfaceQpValueUserObject(const InputParameters & parameters)
  : InterfaceQpValueUserObjectBase(parameters),
    _u(_use_old_value ? coupledValueOld("var") : coupledValue("var")),
    _u_neighbor(
        parameters.isParamSetByUser("var_neighbor")
            ? (_use_old_value ? coupledNeighborValueOld("var_neighbor")
                              : coupledNeighborValue("var_neighbor"))
            : (_use_old_value ? coupledNeighborValueOld("var") : coupledNeighborValue("var")))

{
}

InterfaceQpValueUserObject::~InterfaceQpValueUserObject() {}

void
InterfaceQpValueUserObject::execute()
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
      // compute average value at qp
      vec[qp].resize(1);
      vec[qp][0] = computeInterfaceValueType(_u[qp], _u_neighbor[qp]);
    }
  }
  else
    mooseError("InterfaceQpValueUserObject:: cannot find the required element " +
               std::to_string(_current_elem->id()) + "and side" + std::to_string(_current_side));
}
