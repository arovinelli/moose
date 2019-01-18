//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef CZMTESTLAW_SVM_H
#define CZMTESTLAW_SVM_H

#include "CZMTractionSeparationUOBase.h"

class CZMTestLaw_SVM;

template <> InputParameters validParams<CZMTestLaw_SVM>();

/**
Traction sepration law basic user object
 */
class CZMTestLaw_SVM : public CZMTractionSeparationUOBase {
public:
  CZMTestLaw_SVM(const InputParameters &parameters);

  RealVectorValue computeTractionLocal(unsigned int qp) const override;
  RankTwoTensor
  computeTractionSpatialDerivativeLocal(unsigned int qp) const override;
  std::vector<RealVectorValue>
  computeTractionOtherVarsDerivatives(unsigned int qp) const override;
};

#endif // CZMTESTLAW_SVM_H
