//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterfaceQpMaterialPropertyRealUO.h"
#include "MooseMesh.h"
registerMooseObject("MooseApp", InterfaceQpMaterialPropertyRealUO);

defineLegacyParams(InterfaceQpMaterialPropertyRealUO);

InputParameters
InterfaceQpMaterialPropertyRealUO::validParams()
{
  InputParameters params = InterfaceQpMaterialPropertyUserObjectBase<Real>::validParams();
  params.addClassDescription("Interface Qp value UO working on Real material properties");
  return params;
}

InterfaceQpMaterialPropertyRealUO::InterfaceQpMaterialPropertyRealUO(
    const InputParameters & parameters)
  : InterfaceQpMaterialPropertyUserObjectBase<Real>(parameters)

{
}

Real
InterfaceQpMaterialPropertyRealUO::computeValueMaster(const unsigned int qp)
{
  std::cout << "_prop[qp]" << _prop[qp] << std::endl;
  if (_compute_rate)
  {
    std::cout << "_prop_old[qp]" << (*_prop_old)[qp] << std::endl;
    return (_prop[qp] - (*_prop_old)[qp]) / _dt;
  }
  else
    return _prop[qp];
}

Real
InterfaceQpMaterialPropertyRealUO::computeValueSlave(const unsigned int qp)
{

  std::cout << "_prop_neighbor[qp]" << _prop_neighbor[qp] << std::endl;
  if (_compute_rate)
  {
    std::cout << "_prop_neighbor_old[qp]" << (*_prop_neighbor_old)[qp] << std::endl;
    return (_prop_neighbor[qp] - (*_prop_neighbor_old)[qp]) / _dt;
  }
  else
    return _prop_neighbor[qp];
}
