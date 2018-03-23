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

#include "NormalTractionAux.h"
#include "Assembly.h"

registerMooseObject("TensorMechanicsApp", NormalTractionAux);

template <>
InputParameters
validParams<NormalTractionAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<MaterialPropertyName>("stress",
        "the variable including the cauchy stress tensor");
  // params.addRequiredCoupledVar("disp_x", "the X component of the displamcent");
  // params.addRequiredCoupledVar("disp_y", "the Y component of the displamcent");
  // params.addRequiredCoupledVar("disp_z", "the Z component of the displamcent");

  return params;
}



NormalTractionAux::NormalTractionAux(const InputParameters & parameters)
  : AuxKernel(parameters),
   _normals(_assembly.normals()),
   _stress(getMaterialProperty<RankTwoTensor>("stress"))
   // _ux(coupledValue("disp_x")),
   // _uy(coupledValue("disp_y")),
   // _uz(coupledValue("disp_z")),

{
}

Real
NormalTractionAux::computeValue()
{
  // std::vector<Real> _displ = {_ux[_qp],_uy[_qp],_uz[_qp]};
  _tNormal = 0.0;
  for (unsigned int i = 0; i < 3; i++ ) {
      for (unsigned int j = 0; j < 3; j++ ) {
  _tNormal += _normals[_qp](j)*_stress[_qp](i,j)*_normals[_qp](i);
}
}

  return _tNormal;
}
