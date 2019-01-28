//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CZMLawShamNeedleman.h"
#include "Material.h"
#include "MooseError.h"
#include "libmesh/utility.h"

registerMooseObject("MooseApp", CZMLawShamNeedleman);
template <>
InputParameters
validParams<CZMLawShamNeedleman>()
{
  InputParameters params = validParams<CZMTractionSeparationUOBase>();
  params.addClassDescription("Testing a nonlocal cohseive law model based oon "
                             "mises stress, no damage");

  params.addParam<unsigned int>("n_stateful_mp", 2, "number of stateful material properties");
  params.addParam<std::vector<std::string>>(
      "stateful_mp_names",
      std::vector<std::string>{"a", "N"},
      "name of stateful material properties, namely a(cavity radius) and "
      "N(cavity density). Use this paramter if you wnat to change theri names");
  params.addParam<std::vector<unsigned int>>("stateful_mp_sizes",
                                             std::vector<unsigned int>(2, 1),
                                             "size of each stateful material properties");
  params.addRequiredParam<std::vector<Real>>(
      "stateful_mp_initial_values",
      "intial stateful material proeprties values, namely a0 and NI. No "
      "default value are assumed");
  params.addParam<unsigned int>(
      "n_non_stateful_mp", 2, "number of NON-stateful material properties");
  params.addParam<std::vector<std::string>>("non_stateful_mp_names",
                                            std::vector<std::string>{"b", "C"},
                                            "name of NON stateful material properties");
  params.addParam<std::vector<unsigned int>>("non_stateful_mp_sizes",
                                             std::vector<unsigned int>(0),
                                             "size of each NON stateful material properties");
  params.addRequiredParam<MaterialPropertyName>(
      "interface_avg_mises_stress_mp",
      "name of the material property representing the average mises stress");
  params.addRequiredParam<MaterialPropertyName>(
      "interface_avg_hydrostatic_stress_mp",
      "name of the material property representing the average hydrostatic "
      "stress");
  params.addRequiredParam<MaterialPropertyName>(
      "interface_avg_eq_strain",
      "name of the material property representing the average equivalent "
      "plastic strain");
  params.addRequiredParam<MaterialPropertyName>(
      "interface_avg_eq_strain_rate",
      "name of the material property representing the average equivalent "
      "plastic strain rate");
  params.addRequiredParam<Real>("E_interface", "interface Young modulus");
  params.addRequiredParam<Real>("beta_exponent", "Traction exponent");
  params.addRequiredParam<Real>("D_gb", "gran boundary diffusion coefficient");
  params.addRequiredParam<Real>("n_exp", "power law creep exponent");
  params.addRequiredParam<Real>("FN", "nucleation rate constant");
  params.addRequiredParam<Real>("Nmax", "initial cavity number density");
  params.addRequiredParam<Real>("psi_angle", "equilibrium cavity tip half-angle [degree]");
  params.addRequiredParam<Real>("sigma_0", "traction normalization parameter");
  params.addRequiredParam<Real>("S_thr", "threshold value of S for nucleation on set");
  params.addParam<bool>("use_old_avg_prop",
                        true,
                        "if true (default) old matrial property will "
                        "be used for the average interface quantities");

  return params;
}

CZMLawShamNeedleman::CZMLawShamNeedleman(const InputParameters & parameters)
  : CZMTractionSeparationUOBase(parameters),
    _use_old_avg_prop(getParam<bool>("use_old_avg_prop")),
    _avg_mises_stress(_use_old_avg_prop
                          ? getMaterialPropertyOld<Real>("interface_avg_mises_stress_mp")
                          : getMaterialProperty<Real>("interface_avg_mises_stress_mp")),
    _avg_hyd_stress(_use_old_avg_prop
                        ? getMaterialPropertyOld<Real>("interface_avg_hydrostatic_stress_mp")
                        : getMaterialProperty<Real>("interface_avg_hydrostatic_stress_mp")),
    _avg_eq_strain(_use_old_avg_prop ? getMaterialPropertyOld<Real>("interface_avg_eq_strain")
                                     : getMaterialProperty<Real>("interface_avg_eq_strain")),
    _avg_eq_strain_rate(_use_old_avg_prop
                            ? getMaterialPropertyOld<Real>("interface_avg_eq_strain_rate")
                            : getMaterialProperty<Real>("interface_avg_eq_strain_rate")),
    _a_old(getMaterialPropertyOldByName<std::vector<Real>>(_stateful_mp_names[0])),
    _N_old(getMaterialPropertyOldByName<std::vector<Real>>(_stateful_mp_names[1])),
    _pi(libMesh::pi),
    _E_interface(getParam<Real>("E_interface")),
    _beta_exponent(getParam<Real>("beta_exponent")),
    _D_gb(getParam<Real>("D_gb")),
    _n_exp(getParam<Real>("n_exp")),
    _psi_angle(getParam<Real>("psi_angle") * _pi / 180.),
    _h((1. / (1. + std::cos(_psi_angle)) - std::cos(_psi_angle) / 2.) / std::sin(_psi_angle)),
    _sigma_0(getParam<Real>("sigma_0")),
    _S_thr(getParam<Real>("S_thr")),
    _FN(getParam<Real>("FN")),
    _NI(getParam<std::vector<Real>>("stateful_mp_initial_values")[1]),
    _Nmax(getParam<Real>("Nmax")),
    _a0(getParam<std::vector<Real>>("stateful_mp_initial_values")[0]),
    _b0(fun_b(_NI)),
    _b_sat(fun_b(_Nmax)),
    _n_vars_newton(3),
    _J_newton(_n_vars_newton, _n_vars_newton),
    _R_newton(_n_vars_newton),
    _X_newton(_n_vars_newton)

{
}

RealVectorValue
CZMLawShamNeedleman::computeTractionLocal(unsigned int qp) const
{
  RealVectorValue TractionLocal;
  // Real avg_svm = (*_other_scalar_mp[0])[qp];
  // TractionLocal(0) = avg_svm * _displacement_jump[qp](0);
  // TractionLocal(1) = 0;
  // TractionLocal(2) = 0;

  return TractionLocal;
}

RankTwoTensor
CZMLawShamNeedleman::computeTractionSpatialDerivativeLocal(unsigned int qp) const
{
  RankTwoTensor TractionSpatialDerivativeLocal;
  // Real avg_svm = (*_other_scalar_mp[0])[qp];
  // for (unsigned int i = 0; i < 3; i++)
  //   for (unsigned int j = 0; j < 3; j++)
  //     TractionSpatialDerivativeLocal(i, j) = 0;
  //
  // TractionSpatialDerivativeLocal(0, 0) = avg_svm;

  return TractionSpatialDerivativeLocal;
}

// std::vector<RealVectorValue>
// CZMLawShamNeedleman::computeTractionOtherVarsDerivatives(
//     unsigned int qp) const {
//   std::vector<RealVectorValue> TractionOtherVarsDerivatives(1,
//                                                             RealVectorValue(0));
//   // for (unsigned int i = 0; i < 1; i++)
//   //   for (unsigned int j = 0; j < 3; j++)
//   //     TractionOtherVarsDerivatives[i](j) = 0;
//   //
//   // TractionOtherVarsDerivatives[0](0) = _displacement_jump[qp](0);
//
//   return TractionOtherVarsDerivatives;
// }

void
CZMLawShamNeedleman::NewtonStep(const std::vector<Real> x_old,
                                std::vector<Real> x_new,
                                std::vector<Real> residaul,
                                const std::vector<std::vector<Real>> J)
{
  // compute the update value of the internal variables
  for (unsigned int i = 0; i < _n_vars_newton; ++i)
    for (unsigned int j = 0; j < _n_vars_newton; ++j)
      _J_newton(i, j) = J[i][j];

  for (unsigned int i = 0; i < _n_vars_newton; ++i)
    _R_newton(i) = residaul[i];

  _J_newton.lu_solve(_R_newton, _X_newton);

  for (unsigned int i = 0; i < _n_vars_newton; ++i)
    x_new[i] = x_old[i] - _X_newton(i);
}
