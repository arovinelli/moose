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

#include "TractionAuxField.h"
#include "Assembly.h"

registerMooseObject("TensorMechanicsApp", TractionAuxField);

template <>
InputParameters
validParams<TractionAuxField>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredParam<unsigned int>("tIdx",
     "the variable including the cauchy stress tensor");
  params.addRequiredCoupledVar("s00", "stress 00");
  params.addRequiredCoupledVar("s01", "stress 01");
  params.addRequiredCoupledVar("s02", "stress 02");
  params.addRequiredCoupledVar("s10", "stress 10");
  params.addRequiredCoupledVar("s11", "stress 11");
  params.addRequiredCoupledVar("s12", "stress 12");
  params.addRequiredCoupledVar("s20", "stress 20");
  params.addRequiredCoupledVar("s21", "stress 21");
  params.addRequiredCoupledVar("s22", "stress 22");
  return params;
}



TractionAuxField::TractionAuxField(const InputParameters & parameters)
  : AuxKernel(parameters),
   _tIdx(getParam<unsigned int>("tIdx")),
   _normals(_assembly.normals()),
   _s00(coupledValue("s00")),
   _s01(coupledValue("s01")),
   _s02(coupledValue("s02")),
   _s10(coupledValue("s10")),
   _s11(coupledValue("s11")),
   _s12(coupledValue("s12")),
   _s20(coupledValue("s20")),
   _s21(coupledValue("s21")),
   _s22(coupledValue("s22"))

   // _traction(getMaterialProperty<std::vector<Real>>("traction"))

{
}

Real
TractionAuxField::computeValue()
{
  _tComp  = 0.;
  const std::vector< std::vector<Real>> _tensorTemp {
                        {_s00[_qp],_s01[_qp],_s02[_qp]},
                        {_s10[_qp],_s11[_qp],_s12[_qp]},
                        {_s20[_qp],_s21[_qp],_s22[_qp]} };

  for (unsigned int j = 0; j < 3; j++ ) {
    _tComp += _tensorTemp[_tIdx][j]*_normals[_qp](j);
  }
   return _tComp;
}
