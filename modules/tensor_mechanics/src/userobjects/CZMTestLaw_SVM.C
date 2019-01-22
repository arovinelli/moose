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

registerMooseObject("MooseApp", CZMTestLaw_SVM);
template <>
InputParameters
validParams<CZMTestLaw_SVM>()
{
  InputParameters params = validParams<CZMTractionSeparationUOBase>();
  params.addClassDescription("Testing a nonlocal cohseive law model based oon "
                             "mises stress, no damage");

  params.addParam<unsigned int>("n_stateful_mp", 0, "number of stateful material properties");
  params.addParam<unsigned int>(
      "n_non_stateful_mp", 0, "number of NON-stateful material properties");
  params.addParam<Real>("interface_stiffness", 5e8, "interface_stiffness");

  return params;
}

CZMTestLaw_SVM::CZMTestLaw_SVM(const InputParameters & parameters)
  : CZMTractionSeparationUOBase(parameters),
    _interface_stiffness(getParam<Real>("interface_stiffness"))

{
}

RealVectorValue
CZMTestLaw_SVM::computeTractionLocal(unsigned int qp) const
{

  RealVectorValue TractionLocal;
  Real s_zz_avg = (*_other_avg_scalar_mp[0])[qp];
  TractionLocal(0) =
      _displacement_jump[qp](0) * _interface_stiffness * (1 - std::pow(s_zz_avg / 10., 2));
  TractionLocal(1) =
      _displacement_jump[qp](1) * _interface_stiffness * (1 - std::pow(s_zz_avg / 10., 2));
  TractionLocal(2) =
      _displacement_jump[qp](2) * _interface_stiffness * (1 - std::pow(s_zz_avg / 10., 2));
  return TractionLocal;
}

RankTwoTensor
CZMTestLaw_SVM::computeTractionSpatialDerivativeLocal(unsigned int qp) const
{

  RankTwoTensor TractionSpatialDerivativeLocal;
  Real s_zz_avg = (*_other_avg_scalar_mp[0])[qp];
  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      if (i == j)
        TractionSpatialDerivativeLocal(i, i) =
            _interface_stiffness * (1 - std::pow(s_zz_avg / 10., 2));
      else
        TractionSpatialDerivativeLocal(i, j) = 0;

  return TractionSpatialDerivativeLocal;
}

RealVectorValue
CZMTestLaw_SVM::computeTractionOtherAveragedScalarVarDerivatives(unsigned int qp,
                                                                 unsigned int mp_index) const
{

  RealVectorValue traction_other_avg_scalar_var_derivatives;
  Real s_zz_avg = (*_other_avg_scalar_mp[0])[qp];
  if (mp_index == 0) /* derivative wrt to vm stress*/
    for (unsigned int j = 0; j < 3; j++)
      traction_other_avg_scalar_var_derivatives(j) =
          _displacement_jump[qp](j) * _interface_stiffness * (-2 * s_zz_avg / std::pow(10., 2));

  return traction_other_avg_scalar_var_derivatives;
}
