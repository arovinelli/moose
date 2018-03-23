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

#ifndef TRACTIONAUX_H
#define TRACTIONAUX_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

class TractionAux;

template <>
InputParameters validParams<TractionAux>();

class TractionAux : public AuxKernel
{
public:

  TractionAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
  const unsigned int _tIdx;
  const MooseArray<Point> & _normals;
  const MaterialProperty<RankTwoTensor> & _stress;
  Real _tComp;
};

#endif // TRACTIONAUX_H
