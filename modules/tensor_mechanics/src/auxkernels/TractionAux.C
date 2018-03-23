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

#include "TractionAux.h"
#include "Assembly.h"

registerMooseObject("TensorMechanicsApp", TractionAux);

template <>
InputParameters
validParams<TractionAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<unsigned int>("tIdx",
     "the variable including the cauchy stress tensor");
  params.addRequiredParam<MaterialPropertyName>("stress",
        "the variable including the cauchy stress tensor");
  return params;
}



TractionAux::TractionAux(const InputParameters & parameters)
  : AuxKernel(parameters),
   _tIdx(getParam<unsigned int>("tIdx")),
   _normals(_assembly.normals()),
   _stress(getMaterialProperty<RankTwoTensor>("stress"))

   // _traction(getMaterialProperty<std::vector<Real>>("traction"))

{
}

Real
TractionAux::computeValue()
{
  _tComp  = 0.;
  for (unsigned int j = 0; j < 3; j++ ) {
    _tComp += _stress[_qp](_tIdx,j)*_normals[_qp](j);
  }
  return _tComp;
}
