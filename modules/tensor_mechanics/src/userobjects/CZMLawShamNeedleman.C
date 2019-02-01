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

  params.addParam<unsigned int>("n_stateful_mp", 3, "number of stateful material properties");
  params.addParam<std::vector<std::string>>(
      "stateful_mp_names",
      std::vector<std::string>{"a", "N", "T"},
      "name of stateful material properties, namely a(cavity radius) and "
      "N(cavity density). Use this paramter if you wnat to change theri names");
  params.addParam<std::vector<unsigned int>>("stateful_mp_sizes",
                                             std::vector<unsigned int>{1, 1, 3},
                                             "size of each stateful material properties");
  params.addRequiredParam<std::vector<Real>>(
      "stateful_mp_initial_values",
      "intial stateful material proeprties values, namely a0 and NI. No "
      "default value are assumed");
  params.addParam<unsigned int>(
      "n_non_stateful_mp", 0, "number of NON-stateful material properties");
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
  params.addParam<Real>("E_penalty", 1, "interface Young modulus");
  params.addRequiredParam<Real>("beta_exponent", "Traction exponent");
  params.addRequiredParam<Real>("D_gb", "gran boundary diffusion coefficient");
  params.addRequiredParam<Real>("n_exponent", "power law creep exponent");
  params.addRequiredParam<Real>("FN", "nucleation rate constant");
  params.addRequiredParam<Real>("Nmax", "initial cavity number density");
  params.addRequiredParam<Real>("psi_angle", "equilibrium cavity tip half-angle [degree]");
  params.addRequiredParam<Real>("sigma_0", "traction normalization parameter");
  params.addRequiredParam<Real>("S_thr", "threshold value of S for nucleation on set");
  params.addParam<bool>("use_old_avg_prop",
                        true,
                        "if true (default) old matrial property will "
                        "be used for the average interface quantities");
  params.addParam<bool>("automatic_scale_opening_variable",
                        true,
                        "automatically scale the variables while solving for the opening mode");

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
    _T_old(getMaterialPropertyOldByName<std::vector<Real>>(_stateful_mp_names[2])),
    _pi(libMesh::pi),
    _E_interface(getParam<Real>("E_interface")),
    _E_penalty(getParam<Real>("E_penalty")),
    _beta_exponent(getParam<Real>("beta_exponent")),
    _D_gb(getParam<Real>("D_gb")),
    _n_exponent(getParam<Real>("n_exponent")),
    _alpha_n(3. / (_n_exponent * 2.)),
    _psi_angle(getParam<Real>("psi_angle") * _pi / 180.),
    _h((1. / (1. + std::cos(_psi_angle)) - std::cos(_psi_angle) / 2.) / std::sin(_psi_angle)),
    _sigma_0(getParam<Real>("sigma_0")),
    _S_thr(getParam<Real>("S_thr")),
    _FN(getParam<Real>("FN")),
    _NI(getParam<std::vector<Real>>("stateful_mp_initial_values")[1]),
    _Nmax(getParam<Real>("Nmax")),
    _a0(getParam<std::vector<Real>>("stateful_mp_initial_values")[0]),
    _b0(fun_b_init(_NI)),
    _b_sat(fun_b_init(_Nmax)),
    // initialize varaibles related to solving for the opening mode

    _var_names_newton_opening({"a", "N", "TN"}),
    _n_vars_newton_opening(_var_names_newton_opening.size()),
    _other_var_names_newton_opening({"b", "S", "duN", "dVi", "dec", "h", "sVM", "sH"}),
    _n_other_var_names_newton_opening(_other_var_names_newton_opening.size()),
    _n_var_names_newton_opening_total(_n_vars_newton_opening + _n_other_var_names_newton_opening),
    _all_var_names_newton_opening(_n_var_names_newton_opening_total),
    _automatic_scale_opening_variable(getParam<bool>("automatic_scale_opening_variable")),
    _standard_scale_factor(_n_vars_newton_opening, 1.),
    _opening_mode_vars(_var_names_newton_opening, _other_var_names_newton_opening, 0, 0)
{
}
//                    std::vector<Real>(3, 0),
//                    _other_var_names_newton_opening,
//                    _standard_scale_factor)

// _n_var_names_newton_opening(_var_names_newton_opening.size()),
// _newton_var_opening(_var_names_newton_opening),
// _R_newton_opening(_n_var_names_newton_opening),
// _J_newton_opening(_n_var_names_newton_opening, _n_var_names_newton_opening),

// void
// CZMLawShamNeedleman::initNLVar()
// {
//   DenseVector<Real> opening_var_initial_values(_n_vars_newton_opening);
//   opening_var_initial_values(0) = 0;
//   opening_var_initial_values(1) = 0;
//   opening_var_initial_values(2) = 0;
//   std::vector<Real> _standard_scale_factor(_n_vars_newton_opening, 1.);
//   _opening_mode_vars(_var_names_newton_opening,
//                      opening_var_initial_values,
//                      _other_var_names_newton_opening,
//                      _standard_scale_factor);
// }

RealVectorValue
CZMLawShamNeedleman::computeTractionLocal(unsigned int qp) const
{
  RealVectorValue TractionLocal;
  // DenseVector<Real> x_old(_n_vars_newton_opening);
  // DenseVector<Real> x_new(_n_vars_newton_opening);
  // DenseVector<Real> residual(_n_vars_newton_opening);
  // DenseMatrix<Real> J(_n_vars_newton_opening, _n_vars_newton_opening);
  // // initialize newton first guess to previous variables values
  // x_old(0) = _a_old[qp][0];
  // x_old(1) = _N_old[qp][0];
  // x_old(2) = _T_old[qp][0]; /* normal traction component*/
  // _displacement_jump[qp](0);
  //
  // /*first newton iteration*/
  // computeNewtonResidualAndJacobian(x_old, residual, J);
  // NewtonStep(x_old, x_new, residual, J);

  // for (unsigned int i = 0; i < _n_vars_newton_opening; ++i)
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
CZMLawShamNeedleman::computeNewtonResidualAndJacobian()
{
  fun_b_N();
  TnElastic_Res(); // residual and Jacobian for Normal Tractyion equation

  DenseVector<Real> R(_n_vars_newton_opening);
  DenseMatrix<Real> J(_n_vars_newton_opening, _n_vars_newton_opening);

  _opening_mode_vars.setNLSResidual(_var_names_newton_opening, R);
  _opening_mode_vars.setNLSJacobian(_var_names_newton_opening, J);
}

void
CZMLawShamNeedleman::Newton()
{
  // check if scaling need to be used to solve the NL system

  // DenseVector<Real> x_old = x_init;
  // DenseVector<Real> x_new = x_init;
  // DenseMatrix<Real> J(_n_vars_newton_opening, _n_vars_newton_opening);
  // DenseVector<Real> R(_n_vars_newton_opening);
  //
  // unsigned int N = 0;
  // Real current_max_error = 1.;
  // newton_flag = 1;
  // do
  // {
  //   x_old = x_new;
  //   computeNewtonResidualAndJacobian(x_old, R, J);
  //   // compute current error
  //
  //   // compute step
  //   NewtonStep(x_old, x_new, R, J);
  //   N++;
  // } while (N < _max_newton_iter && current_max_error > _max_abs_tol);
  // if (current_max_error < _max_abs_tol)
  // {
  //   newton_flag = 0;
  //   x_sol = x_old;
  // need to do somehitng about the dx_dDeltau
}

void
CZMLawShamNeedleman::NewtonStep()
{
  DenseMatrix<Real> J_temp = _opening_mode_vars.getNLSJacobianScaled();
  DenseVector<Real> R_temp = _opening_mode_vars.getNLSResidual();
  DenseVector<Real> x_old = _opening_mode_vars.getNLSValueScaledOld();
  DenseVector<Real> x_new(_n_vars_newton_opening);

  // compute the update value of the internal variables
  J_temp.lu_solve(R_temp, x_new);

  for (unsigned int i = 0; i < _n_vars_newton_opening; ++i)
    x_new(i) =
        (x_old(i) - x_new(i)) * _opening_mode_vars.getScaleFactor(_var_names_newton_opening[i]);

  _opening_mode_vars.advanceNLS();
  _opening_mode_vars.setNLSValue(_var_names_newton_opening, x_new);
}

void
CZMLawShamNeedleman::computeCompliance()
{
  // Real a = 0, b = 0, N = 0, TN = 0;
  // unrollVars(x_current, a, b, N, TN);
  Real a = _opening_mode_vars.getValue("a");
  Real b = _opening_mode_vars.getValue("b");
  Real uN = getuN();
  Real K = 6;
  Real db_dN = _opening_mode_vars.getDerivative("b", "N");
  //
  Real S = K * _b_sat / (_E_interface * (1. - (a / b)));

  Real dS_da = K * b * _b_sat / (_E_interface * pow(a - b, 2.));
  Real dS_dN = -K * a * _b_sat / (_E_interface * pow(a - b, 2.)) * db_dN;
  //
  if (uN < 0.)
  {
    S /= _E_penalty;
    dS_da /= _E_penalty;
    dS_dN /= _E_penalty;
  }

  _opening_mode_vars.setValue("S", S);
  _opening_mode_vars.setDerivative("S", "a", dS_da);
  _opening_mode_vars.setDerivative("S", "N", dS_dN);
}

void
CZMLawShamNeedleman::Tn_fun()
{
  computeCompliance();
  Real N = _opening_mode_vars.getValue("N");
  Real Tn_prev = _opening_mode_vars.getInitialValue("TN");
  Real duN = getDeltauN(); // between the current and previous solution
  Real dVi = _opening_mode_vars.getValue("dVi");
  DenseVector<Real> ddVi_dx = _opening_mode_vars.getDerivativeVector("dVi");

  Real S = _opening_mode_vars.getValue("S");
  DenseVector<Real> dS_dx = _opening_mode_vars.getDerivativeVector("S");

  Real Num = duN - dVi * N;
  Real Tn_comp = Num / S + Tn_prev;

  _opening_mode_vars.setComputedValue("TN", Tn_comp);

  DenseVector<Real> dN_dx = _opening_mode_vars.getUnitDerivativeVector("N");
  DenseVector<Real> dduN_dx = _opening_mode_vars.getUnitDerivativeVector("duN");

  for (unsigned int i = 0; i < _n_var_names_newton_opening_total; i++)
  {
    std::string var_name = _all_var_names_newton_opening[i];
    Real dNum_dx = dduN_dx(i) - (ddVi_dx(i) * N + dN_dx(i) * dVi);
    _opening_mode_vars.setDerivative(
        "TN", var_name, (dNum_dx * S - Num * dS_dx(i)) / std::pow(S, 2.));
  }
}

void
CZMLawShamNeedleman::TnElastic_Res()
{
  Tn_fun();
  Real R_TN = _opening_mode_vars.getValue("TN") - _opening_mode_vars.getComputedValue("TN");
  _opening_mode_vars.setResidual("TN", R_TN);
  DenseVector<Real> dTnComp_dx = _opening_mode_vars.getDerivativeVector("TN");
  DenseVector<Real> dTn_dx = _opening_mode_vars.getUnitDerivativeVector("TN");

  for (unsigned int i = 0; i < _n_var_names_newton_opening_total; i++)
  {
    std::string var_name = _all_var_names_newton_opening[i];
    _opening_mode_vars.setDerivative("TN", var_name, dTn_dx(i) + dTnComp_dx(i));
  }
}

void
CZMLawShamNeedleman::fun_b_N()
{
  Real N = _opening_mode_vars.getValue("N");
  Real b = fun_b_init(N);
  Real db_dN = -0.5 * std::pow(_pi * N, -1.5) * _pi;
  _opening_mode_vars.setValue("b", b);
  _opening_mode_vars.setDerivative("b", "N", db_dN);
}

void
CZMLawShamNeedleman::N_Fun()
{
  unsigned int qp = _opening_mode_vars.get_qp();
  Real Tn = _opening_mode_vars.getValue("TN");
  Real N_prev = _opening_mode_vars.getInitialValue("N");
  Real N = _opening_mode_vars.getValue("N");
  Real fVal = N_prev;
  Real dN_dTN = 0;
  Real dN_dec = 0;

  if (Tn >= 0 && N <= _Nmax && N >= N_prev)
  {

    Real deltaEC = getDeltaEc();
    Real S = std::pow((Tn / _sigma_0), _beta_exponent) * _avg_eq_strain[qp];
    if (S > _S_thr && deltaEC > 0)
    {
      Real delatN = _FN * std::pow(Tn / _sigma_0, _beta_exponent) * deltaEC;
      Real N = delatN + N_prev;
      if (N > N_prev)
      {
        fVal = N;
        dN_dec = _FN * std::pow(Tn / _sigma_0, _beta_exponent);
        dN_dTN = _beta_exponent * _FN * deltaEC * std::pow(Tn / _sigma_0, _beta_exponent - 1.) /
                 _sigma_0;
      }
    }
  }

  _opening_mode_vars.setComputedValue("N", fVal);
  _opening_mode_vars.setDerivative("N", "TN", dN_dTN);
  _opening_mode_vars.setDerivative("N", "dec", dN_dec);
}

void
CZMLawShamNeedleman::N_Res()
{
  N_Fun();
  Real N = _opening_mode_vars.getValue("N");
  Real NComp = _opening_mode_vars.getComputedValue("N");
  DenseVector<Real> dNComp_dx = _opening_mode_vars.getDerivativeVector("N");
  DenseVector<Real> dN_dx = _opening_mode_vars.getUnitDerivativeVector("N");

  Real R_N = N - NComp;
  _opening_mode_vars.setResidual("N", R_N);
  for (unsigned int i = 0; i < _n_var_names_newton_opening_total; i++)
  {
    std::string var_name = _all_var_names_newton_opening[i];
    _opening_mode_vars.setDerivative("N", var_name, dN_dx(i) - dNComp_dx(i));
  }
}

void
CZMLawShamNeedleman::a_Fun()
{
  Real dV = _opening_mode_vars.getValue("dVi");
  DenseVector<Real> ddV_dx = _opening_mode_vars.getDerivativeVector("dVi");
  Real a_prev = _opening_mode_vars.getInitialValue("a");

  Real a = _opening_mode_vars.getValue("a");

  Real denominator = 4. * _pi * _h * std::pow(a, 2.);
  DenseVector<Real> dden_dx = _opening_mode_vars.getZeroDerivativeVector();

  setZeroVectorDerivative(dden_dx, "a", 8 * _pi * _h * a);
  setZeroVectorDerivative(dden_dx, "h", 4 * _pi * std::pow(a, 2.));

  Real aComp = dV / denominator + a_prev;
  _opening_mode_vars.setComputedValue("a", aComp);

  for (unsigned int i = 0; i < _n_var_names_newton_opening_total; i++)
  {
    std::string var_name = _all_var_names_newton_opening[i];
    Real da_dx = (ddV_dx(i) * denominator - dV * dden_dx(i)) / std::pow(denominator, 2.);
    _opening_mode_vars.setDerivative("a", var_name, da_dx);
  }
}

void
CZMLawShamNeedleman::a_Res()
{
  a_Fun();
  Real a = _opening_mode_vars.getValue("a");
  Real aComp = _opening_mode_vars.getComputedValue("a");
  DenseVector<Real> daComp_dx = _opening_mode_vars.getDerivativeVector("a");
  Real R_a = a - aComp;
  _opening_mode_vars.setResidual("a", R_a);
  DenseVector<Real> da_dx = _opening_mode_vars.getUnitDerivativeVector("a");

  for (unsigned int i = 0; i < _n_var_names_newton_opening_total; i++)
  {
    std::string var_name = _all_var_names_newton_opening[i];
    _opening_mode_vars.setDerivative("a", var_name, da_dx(i) - daComp_dx(i));
  }
}

void
CZMLawShamNeedleman::dVi_Fun()
{
  Real dvl1, dvl2, dvh1, dvh2, dVL, dVH, dVi;
  DenseVector<Real> ddvl1_dx = _opening_mode_vars.getZeroDerivativeVector();
  DenseVector<Real> ddvl2_dx = _opening_mode_vars.getZeroDerivativeVector();
  DenseVector<Real> ddvh1_dx = _opening_mode_vars.getZeroDerivativeVector();
  DenseVector<Real> ddvh2_dx = _opening_mode_vars.getZeroDerivativeVector();
  DenseVector<Real> ddVi_dx = _opening_mode_vars.getZeroDerivativeVector();

  dV_1_Fun(dvl1, ddvl1_dx, 1);
  dV_2_Fun(dvl2, ddvl2_dx, 1);
  dV_1_Fun(dvh1, ddvh1_dx, 2);
  dV_2_Fun(dvh2, ddvh2_dx, 2);

  dVL = dvl1 + dvl2;
  dVH = dvh1 + dvh2;

  if (std::abs(dVL) >= std::abs(dVH))
  {
    dVi = dVL;
    for (unsigned int i = 0; i < _n_var_names_newton_opening_total; i++)
      ddVi_dx(i) = ddvl1_dx(i) + ddvl2_dx(i);
  }
  else if (std::abs(dVH) > std::abs(dVL))
  {
    dVi = dVH;
    for (unsigned int i = 0; i < _n_var_names_newton_opening_total; i++)
      ddVi_dx(i) = ddvh1_dx(i) + ddvh2_dx(i);
  }
  else
    mooseError(
        "dVi_Fun: cannot determine the max between dVL1 adn dVH1: dVL1 = ", dVL, " dVH1 ", dVH);

  _opening_mode_vars.setValue("dVi", dVi);

  for (unsigned int i = 0; i < _n_var_names_newton_opening_total; i++)
  {
    std::string var_name = _all_var_names_newton_opening[i];
    _opening_mode_vars.setDerivative("dVi", var_name, ddVi_dx(i));
  }
}

void
CZMLawShamNeedleman::dV_1_Fun(Real & dv1, DenseVector<Real> & ddv1_dx, const unsigned int i_vol)
{
  Real uN = getuN();
  Real a = _opening_mode_vars.getValue("a");
  Real Tn = _opening_mode_vars.getValue("TN");
  Real dt = _opening_mode_vars.getDt();
  if (uN >= 0 && ((a > _a0 && Tn < 0) || Tn >= 0))
  {
    Real q;
    DenseVector<Real> dq_dx = _opening_mode_vars.getZeroDerivativeVector();

    qi_Fun(q, dq_dx, i_vol);

    Real C = 8. * _pi * _D_gb * dt;
    dv1 = C * Tn / q;
    Real ddv1_dq = -C * Tn / std::pow(q, 2.);

    setZeroVectorDerivative(ddv1_dx, "TN", C / q);
    for (unsigned int i = 0; i < _n_var_names_newton_opening_total; i++)
    {
      std::string var_name = _all_var_names_newton_opening[i];
      setZeroVectorDerivative(ddv1_dx, var_name, ddv1_dx(i) + ddv1_dq * dq_dx(i));
    }
  }
  else
    dv1 = 0.;
}

void
CZMLawShamNeedleman::qi_Fun(Real & q, DenseVector<Real> & dq_dx, const unsigned int i_vol)
{

  Real f = 0;
  DenseVector<Real> df_dx = _opening_mode_vars.getZeroDerivativeVector();

  fi_Fun(f, df_dx, i_vol);
  q = 2. * std::log(1. / f) - (1. - f) * (3. - f);
  Real dq_df = -2. / f - 2. * f + 4.;

  for (unsigned int i = 0; i < _n_var_names_newton_opening_total; i++)
  {
    std::string var_name = _all_var_names_newton_opening[i];
    setZeroVectorDerivative(dq_dx, var_name, dq_df * df_dx(i));
  }
}

void
CZMLawShamNeedleman::fL_Fun(Real & fL, DenseVector<Real> & dfL_dx)
{
  Real f_ab;
  DenseVector<Real> f_ab_dx = _opening_mode_vars.getZeroDerivativeVector();
  Real f_abL;
  DenseVector<Real> f_abL_dx = _opening_mode_vars.getZeroDerivativeVector();
  f_a_over_b_square_Fun(f_ab, f_ab_dx);
  f_abL_Fun(f_abL, f_abL_dx);

  if (f_ab >= f_abL)
  {
    fL = f_ab;
    dfL_dx = f_ab_dx;
  }
  else if (f_abL > f_ab)
  {
    fL = f_abL;
    dfL_dx = f_abL_dx;
  }
  else
    mooseError(
        "fL_Fun: cannot determine the max between f_ab adn f_abL: f_ab = ", f_ab, " f_abL ", f_abL);
}

void
CZMLawShamNeedleman::fH_Fun(Real & fH, DenseVector<Real> & dfH_dx)
{
  f_a_over_b_square_Fun(fH, dfH_dx);
}

void
CZMLawShamNeedleman::f_a_over_b_square_Fun(Real & f_ab, DenseVector<Real> & df_ab_dx)
{
  Real a = _opening_mode_vars.getValue("a");
  Real b = _opening_mode_vars.getValue("b");
  DenseVector<Real> db_dx = _opening_mode_vars.getDerivativeVector("b");

  f_ab = std::pow(a / b, 2.);
  Real df_ab_da = 2. * a / std::pow(b, 2.);
  setZeroVectorDerivative(df_ab_dx, "a", df_ab_da);

  Real df_ab_db = -2. * std::pow(a, 2.) / std::pow(b, 3.);
  for (unsigned int i = 0; i < _n_var_names_newton_opening_total; i++)
  {
    std::string var_name = _all_var_names_newton_opening[i];
    setZeroVectorDerivative(df_ab_dx, var_name, df_ab_dx(i) + df_ab_db * db_dx(i));
  }
}

void
CZMLawShamNeedleman::f_abL_Fun(Real & f_abL, DenseVector<Real> & df_abL_dx)
{

  Real dec = getDeltaEc();
  f_abL = 0.;
  Real sign_dec = std::copysign(1.0, dec);
  dec = std::abs(dec);

  if (dec > 0.)
  {
    Real a = _opening_mode_vars.getValue("a");
    Real dt = _opening_mode_vars.getDt();
    Real sVM = getAvgSVM();
    DenseVector<Real> dL_dx = _opening_mode_vars.getZeroDerivativeVector();

    Real n = 1. / 3.;
    Real L = sign_dec * std::pow(_D_gb * sVM * dt / (dec), n);
    Real dL_dsVM = sign_dec * n * L / sVM;
    Real dL_ddec = -n * L / (dec);

    setZeroVectorDerivative(dL_dx, "sVM", dL_dsVM);
    setZeroVectorDerivative(dL_dx, "dec", dL_ddec);

    Real a2 = std::pow(a, 2.);
    f_abL = a2 / std::pow(a + 1.5 * L, 2.);

    Real df_abL_da = 3. * (L * a) / std::pow(1.5 * L + a, 3.);

    setZeroVectorDerivative(df_abL_dx, "a", df_abL_da);
    Real df_abL_dL = -3. * a2 / std::pow(1.5 * L + a, 3.);

    for (unsigned int i = 0; i < _n_var_names_newton_opening_total; i++)
    {
      std::string var_name = _all_var_names_newton_opening[i];
      setZeroVectorDerivative(df_abL_dx, var_name, df_abL_dx(i) + df_abL_dL * dL_dx(i));
    }
  }
}

void
CZMLawShamNeedleman::fi_Fun(Real & f, DenseVector<Real> & df_dx, const unsigned int i_vol)
{
  if (i_vol == 1)
    fL_Fun(f, df_dx);
  else if (i_vol == 2)
    fH_Fun(f, df_dx);
  else
    mooseError("Fi_Fun: the supplied value of i_vol is wrong, ", i_vol);
}

void
CZMLawShamNeedleman::dV_2_Fun(Real & dv2, DenseVector<Real> & ddv2_dx, const unsigned int i_vol)
{

  dv2 = 0; // if traction is negative non volume change
  if (_opening_mode_vars.getValue("TN") >= 0)
  {
    if (i_vol == 1)
      dVL2_Fun(dv2, ddv2_dx);
    else if (i_vol == 2)
      dVH2_Fun(dv2, ddv2_dx);
    else
      mooseError("dV_2_Fun: wrong ivol_paramter: i_vol ", i_vol);
  }
}

void
CZMLawShamNeedleman::dVL2_Fun(Real & dvl2, DenseVector<Real> & ddvl2_dx)
{
  dvl2 = 0.;
  if (_opening_mode_vars.getValue("TN") >= 0)
  {
    Real sH = getAvgSH();
    Real sVM = getAvgSVM();
    Real dec = getDeltaEc();

    Real a = _opening_mode_vars.getValue("a");

    Real n = _n_exponent;
    Real h = _h;

    Real m = m_Fun();
    Real beta_n_m = beta_n_m_Fun();
    Real abs_sH_sVM = std::abs(sH / sVM);

    Real ddvl2_da, ddvl2_dec, ddvl2_dSH, ddvl2_dSVM;

    if (abs_sH_sVM >= 1.)
    {
      dvl2 = 2. * _pi * std::pow(a, 3.) * dec * h * m *
             std::pow(_alpha_n * std::abs(sH) / sVM + beta_n_m, n);
      ddvl2_da = 6. * _pi * std::pow(a, 2.) * dec * h * m *
                 std::pow(_alpha_n * std::abs(sH) / sVM + beta_n_m, n);
      ddvl2_dec = 2. * _pi * std::pow(a, 3.) * h * m *
                  std::pow(_alpha_n * std::abs(sH) / sVM + beta_n_m, n);
      ddvl2_dSH = 2. * _pi * std::pow(a, 3.) * _alpha_n * dec * h * m * n * std::pow(sVM, -n) *
                  std::pow(_alpha_n * std::abs(sH) + beta_n_m * sVM, n - 1.) * m;
      ddvl2_dSVM = -2. * _pi * std::pow(a, 3.) * _alpha_n * dec * h * m * n *
                   std::pow(sVM, -n - 1.) *
                   std::pow(_alpha_n * std::abs(sH) + beta_n_m * sVM, n - 1.) * std::abs(sH);
    }
    else
    {
      dvl2 = 2. * _pi * std::pow(a, 3.) * dec * h * sH * std::pow(_alpha_n + beta_n_m, n) / sVM;
      ddvl2_da = 6. * _pi * std::pow(a, 2.) * dec * h * sH * std::pow(_alpha_n + beta_n_m, n) / sVM;
      ddvl2_dec = 2. * _pi * std::pow(a, 3.) * h * sH * std::pow(_alpha_n + beta_n_m, n) / sVM;
      ddvl2_dSH = 2. * _pi * std::pow(a, 3.) * dec * h * std::pow(_alpha_n + beta_n_m, n) / sVM;
      ddvl2_dSVM = -2. * _pi * std::pow(a, 3.) * dec * h * sH * std::pow(_alpha_n + beta_n_m, n) /
                   std::pow(sVM, 2.);
    }

    setZeroVectorDerivative(ddvl2_dx, "a", ddvl2_da);
    setZeroVectorDerivative(ddvl2_dx, "dec", ddvl2_dec);
    setZeroVectorDerivative(ddvl2_dx, "sH", ddvl2_dSH);
    setZeroVectorDerivative(ddvl2_dx, "sVM", ddvl2_dSVM);
  }
  else
    mooseError("dVL2_Fun should not be called because TN < 0: ", _opening_mode_vars.getValue("TN"));
}

void
CZMLawShamNeedleman::dVH2_Fun(Real & dvh2, DenseVector<Real> & ddvh2_dx)
{
  if (_opening_mode_vars.getValue("TN") >= 0)
  {
    Real sH = getAvgSH();
    Real sVM = getAvgSVM();
    Real dec = getDeltaEc();

    Real a = _opening_mode_vars.getValue("a");
    Real b = _opening_mode_vars.getValue("b");
    Real db_dN = _opening_mode_vars.getDerivative("b", "N");
    Real n = _n_exponent;
    Real h = _h;

    Real m = m_Fun();
    // Real beta_n_m = beta_n_m_Fun();
    Real abs_sH_sVM = std::abs(sH / sVM);

    Real ddvh2_da, ddvh2_dN, ddvh2_dec, ddvh2_dSH, ddvh2_dSVM;

    if (abs_sH_sVM >= 1.)
    {
      dvh2 =
          2. * _pi * std::pow(a, 3.) * dec * h * m *
          std::pow((_alpha_n * std::abs(sH) / sVM + m / n) / (-std::pow(0.87 * a / b, 3. / n) + 1.),
                   n);
      ddvh2_da = 6. * _pi * std::pow(a, 2.) * std::pow(b, 3. / n) * dec * h * m *
                 std::pow(std::pow(b, 3. / n) * (_alpha_n * n * std::abs(sH) + m * sVM) /
                              (n * sVM * (pow(b, 3. / n) - pow(0.87 * a, 3. / n))),
                          n) /
                 (std::pow(b, 3. / n) - pow(0.87 * a, 3. / n));
      ddvh2_dN = 6. * _pi * std::pow(a, 3.) * dec * h * m * std::pow(0.87 * a, 3. / n) *
                 std::pow(-std::pow(b, 3. / n) * (_alpha_n * n * std::abs(sH) + m * sVM) /
                              (n * sVM * (-std::pow(b, 3. / n) + std::pow(0.87 * a, 3. / n))),
                          n) /
                 (b * (-std::pow(b, 3. / n) + std::pow(0.87 * a, 3. / n))) * db_dN;
      ddvh2_dec = 2. * _pi * pow(a, 3.) * h * m *
                  std::pow(-std::pow(b, 3. / n) * (_alpha_n * n * std::abs(sH) + m * sVM) /
                               (n * sVM * (-std::pow(b, 3. / n) + std::pow(0.87 * a, 3. / n))),
                           n);
      ddvh2_dSH = 2. * _pi * std::pow(a, 3.) * _alpha_n * dec * h * m * std::pow(n, 2.) *
                  std::pow(-std::pow(b, 3. / n) * (_alpha_n * n * std::abs(sH) + m * sVM) /
                               (n * sVM * (-std::pow(b, 3. / n) + std::pow(0.87 * a, 3. / n))),
                           n) *
                  std::copysign(1, sH) / (_alpha_n * n * std::abs(sH) + m * sVM);
      ddvh2_dSVM = -2. * _pi * std::pow(a, 3.) * _alpha_n * dec * h * m * std::pow(n, 2.) *
                   std::pow(std::pow(b, 3. / n) * (_alpha_n * n * std::abs(sH) + m * sVM) /
                                (n * sVM * (std::pow(b, 3. / n) - pow(0.87 * a, 3. / n))),
                            n) *
                   std::abs(sH) / (sVM * (_alpha_n * n * std::abs(sH) + m * sVM));
    }
    else
    {

      dvh2 =
          2. * _pi * std::pow(a, 3.) * dec * h * m *
          std::pow((_alpha_n * std::abs(sH) / sVM + m / n) / (-std::pow(0.87 * a / b, 3. / n) + 1.),
                   n);
      ddvh2_da = 6. * _pi * std::pow(a, 2.) * std::pow(b, 3. / n) * dec * h * m * sH *
                 std::pow(std::pow(b, 3. / n) * (_alpha_n * n + m) /
                              (n * (std::pow(b, 3. / n) - std::pow(0.87 * a, 3. / n))),
                          n) /
                 (sVM * (std::pow(b, 3. / n) - std::pow(0.87 * a, 3. / n)));
      ddvh2_dN = 6. * _pi * std::pow(a, 3.) * dec * h * m * sH * std::pow(0.87 * a, 3. / n) *
                 std::pow(-std::pow(b, 3. / n) * (_alpha_n * n + m) /
                              (n * (-std::pow(b, 3. / n) + pow(0.87 * a, 3. / n))),
                          n) /
                 (b * sVM * (-std::pow(b, 3. / n) + std::pow(0.87 * a, 3. / n))) * db_dN;
      ddvh2_dec = 2. * _pi * std::pow(a, 3.) * h * m * sH *
                  std::pow(-std::pow(b, 3. / n) * (_alpha_n * n + m) /
                               (n * (-std::pow(b, 3. / n) + std::pow(0.87 * a, 3. / n))),
                           n) /
                  sVM;
      ddvh2_dSH = 2. * _pi * std::pow(a, 3.) * dec * h * m *
                  std::pow(-std::pow(b, 3. / n) * (_alpha_n * n + m) /
                               (n * (-std::pow(b, 3. / n) + std::pow(0.87 * a, 3. / n))),
                           n) /
                  sVM;
      ddvh2_dSVM = -2. * _pi * std::pow(a, 3.) * dec * h * m * sH *
                   std::pow(std::pow(b, 3. / n) * (_alpha_n * n + m) /
                                (n * (std::pow(b, 3. / n) - std::pow(0.87 * a, 3. / n))),
                            n) /
                   std::pow(sVM, 2.);
    }
    setZeroVectorDerivative(ddvh2_dx, "a", ddvh2_da);
    setZeroVectorDerivative(ddvh2_dx, "N", ddvh2_dN);
    setZeroVectorDerivative(ddvh2_dx, "dec", ddvh2_dec);
    setZeroVectorDerivative(ddvh2_dx, "sH", ddvh2_dSH);
    setZeroVectorDerivative(ddvh2_dx, "sVM", ddvh2_dSVM);
  }
  else
    mooseError("dVH2_Fun should not be called because TN < 0: ", _opening_mode_vars.getValue("TN"));
}

int
CZMLawShamNeedleman::m_Fun()
{
  Real sH = getAvgSH();
  return std::copysign(1, sH);
}

Real
CZMLawShamNeedleman::g_Fun()
{
  Real g;
  int m = m_Fun();
  if (m == 1)
    g = std::log(3.) - 2. / 3.;
  else if (m == -1)
    g = 2. * _pi / (9. * std::pow(3., 0.5));
  else
    g = 0.;
  return g;
}

Real
CZMLawShamNeedleman::beta_n_m_Fun()
{
  Real g = g_Fun();
  return (_n_exponent - 1.) * (_n_exponent + g) / pow(_n_exponent, 2.);
}

Real
CZMLawShamNeedleman::ab_ratio(const bool & old /* = false*/)
{
  Real r;
  if (old)
    r = _opening_mode_vars.getValueOld("a") / _opening_mode_vars.getValueOld("b");
  else
    r = _opening_mode_vars.getValue("a") / _opening_mode_vars.getValue("b");
  return r;
}

Real
CZMLawShamNeedleman::getDeltaEc()
{
  Real dt_newton = _opening_mode_vars.getDt();
  unsigned int qp = _opening_mode_vars.get_qp();
  return _avg_eq_strain_rate[qp] * dt_newton;
}

Real
CZMLawShamNeedleman::getEcRate()
{
  unsigned int qp = _opening_mode_vars.get_qp();
  return _avg_eq_strain_rate[qp];
}

Real
CZMLawShamNeedleman::getAvgSVM()
{
  unsigned int qp = _opening_mode_vars.get_qp();
  return _avg_mises_stress[qp];
}

Real
CZMLawShamNeedleman::getAvgSH()
{
  unsigned int qp = _opening_mode_vars.get_qp();
  return _avg_hyd_stress[qp];
}

Real
CZMLawShamNeedleman::getDeltauN()
{
  Real dt_newton = _opening_mode_vars.getDt();
  unsigned int qp = _opening_mode_vars.get_qp();
  return (_displacement_jump[qp](0) - _displacement_jump_old[qp](0)) / _dt * dt_newton;
}

Real
CZMLawShamNeedleman::getuN()
{
  unsigned int qp = _opening_mode_vars.get_qp();
  return _displacement_jump_old[qp](0) + getDeltauN();
}
