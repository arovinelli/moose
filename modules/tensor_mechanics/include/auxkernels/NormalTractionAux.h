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

#ifndef NORMALTRACTIONAUX_H
#define NORMALTRACTIONAUX_H

#include "AuxKernel.h"
#include "RankTwoTensor.h"

class NormalTractionAux;

template <>
InputParameters validParams<NormalTractionAux>();

class NormalTractionAux : public AuxKernel
{
public:

  NormalTractionAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
  const MooseArray<Point> & _normals;
  const MaterialProperty<RankTwoTensor> & _stress;
  Real _tNormal;
};

#endif // TRACTIONAUX_H
