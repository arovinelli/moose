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

#ifndef TRACTIONAUXFIELD_H
#define TRACTIONAUXFIELD_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

class TractionAuxField;

template <>
InputParameters validParams<TractionAuxField>();

class TractionAuxField : public AuxKernel
{
public:

  TractionAuxField(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
  const unsigned int _tIdx;
  const MooseArray<Point> & _normals;
  Real _tComp;
  const VariableValue & _s00;
  const VariableValue & _s01;
  const VariableValue & _s02;
  const VariableValue & _s10;
  const VariableValue & _s11;
  const VariableValue & _s12;
  const VariableValue & _s20;
  const VariableValue & _s21;
  const VariableValue & _s22;
};

#endif // TRACTIONAUXFIELD_H
