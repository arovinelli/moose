//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CZMTestLaw_SVM.h"
#include "Material.h"
#include "MooseError.h"

registerMooseObject("DeerApp", CZMTestLaw_SVM);
template <> InputParameters validParams<CZMTestLaw_SVM>() {
  InputParameters params = validParams<CZMTractionSeparationUOBase>();
  params.addClassDescription("Testing a nonlocal cohseive law model based oon "
                             "mises stress, no damage");

  params.addParam<unsigned int>("n_stateful_mp", 0,
                                "number of stateful material properties");
  params.addParam<unsigned int>("n_non_stateful_mp", 0,
                                "number of NON-stateful material properties");

  return params;
}

CZMTestLaw_SVM::CZMTestLaw_SVM(const InputParameters &parameters)
    : CZMTractionSeparationUOBase(parameters)

{}

RealVectorValue CZMTestLaw_SVM::computeTractionLocal(unsigned int qp) const {

  RealVectorValue TractionLocal;
  Real avg_svm = (*_other_scalar_mp[0])[qp];
  TractionLocal(0) = avg_svm * _displacement_jump[qp](0);
  TractionLocal(1) = 0;
  TractionLocal(2) = 0;

  return TractionLocal;
}

RankTwoTensor
CZMTestLaw_SVM::computeTractionSpatialDerivativeLocal(unsigned int qp) const {

  RankTwoTensor TractionSpatialDerivativeLocal;
  Real avg_svm = (*_other_scalar_mp[0])[qp];
  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      TractionSpatialDerivativeLocal(i, j) = 0;

  TractionSpatialDerivativeLocal(0, 0) = avg_svm;

  return TractionSpatialDerivativeLocal;
}

std::vector<RealVectorValue>
CZMTestLaw_SVM::computeTractionOtherVarsDerivatives(unsigned int qp) const {

  std::vector<RealVectorValue> TractionOtherVarsDerivatives(1,
                                                            RealVectorValue(0));
  for (unsigned int i = 0; i < 1; i++)
    for (unsigned int j = 0; j < 3; j++)
      TractionOtherVarsDerivatives[i](j) = 0;

  TractionOtherVarsDerivatives[0](0) = _displacement_jump[qp](0);

  return TractionOtherVarsDerivatives;
}
