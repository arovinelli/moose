//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef CZMCOPENETRATIONPENALTY_H
#define CZMCOPENETRATIONPENALTY_H

#include "CZMTractionSeparationUOBase.h"

class CZMCopenetrationPenalty;

template <>
InputParameters validParams<CZMCopenetrationPenalty>();

/**
Traction sepration law basic user object
 */
class CZMCopenetrationPenalty : public CZMTractionSeparationUOBase
{
public:
  CZMCopenetrationPenalty(const InputParameters & parameters);

  RealVectorValue computeTractionLocal(unsigned int qp) const override;
  RankTwoTensor computeTractionSpatialDerivativeLocal(unsigned int qp) const override;

  // std::vector<Real> getNewStatefulMaterialProperty(unsigned int qp,
  //                                                  unsigned int mp_index) const override;
  // std::vector<Real> getNewNonStatefulMaterialProperty(unsigned int qp,
  //                                                     unsigned int mp_index) const override;

  // Real getEffectiveJump(unsigned int /*qp*/) const override;

  // unsigned int checkLoadUnload(unsigned int /*qp*/) const override;

protected:
  // cohesive law parameters
  // const Real _displacement_jump_peak;
  // const Real _traction_peak;
  // const Real _beta;

  // const MaterialProperty<std::vector<Real>> & _effective_jump;
  // const MaterialProperty<std::vector<Real>> & _effective_jump_old;

  const MaterialProperty<std::vector<Real>> & _weighted_displacement_jump;
  const MaterialProperty<std::vector<Real>> & _displacement_jump_weights;
  const Real _copenetration_penalty_stiffness;

  // const MaterialProperty<std::vector<Real>> & _effective_traction;

  // Real getEffectiveTraction(unsigned int /*qp*/) const;
  // Real getEffectiveTractionLinear(unsigned int /*qp*/) const;
  // std::vector<Real> getTractionNonLinear(unsigned int /*qp*/) const;
  // std::vector<Real> getTractionLinear(unsigned int /*qp*/) const;
  // std::vector<std::vector<Real>> getTractionSpatialDerivativeNonLinear(unsigned int /*qp*/)
  // const; std::vector<std::vector<Real>> getTractionSpatialDerivativeLinear(unsigned int /*qp*/)
  // const;
};

#endif // CZMCOPENETRATIONPENALTY_H
