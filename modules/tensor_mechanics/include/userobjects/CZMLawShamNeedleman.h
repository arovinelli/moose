//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef CZMLAWSHAMNEEDLEMAN_H
#define CZMLAWSHAMNEEDLEMAN_H

#include "CZMTractionSeparationUOBase.h"
#include "NLSolverVar.h"

class CZMLawShamNeedleman;

template <>
InputParameters validParams<CZMLawShamNeedleman>();

/**
Traction sepration law basic user object
 */
class CZMLawShamNeedleman : public CZMTractionSeparationUOBase
{
public:
  CZMLawShamNeedleman(const InputParameters & parameters);

  RealVectorValue computeTractionLocal(unsigned int qp) const override;
  RankTwoTensor computeTractionSpatialDerivativeLocal(unsigned int qp) const override;
  // std::vector<RealVectorValue>
  // computeTractionOtherVarsDerivatives(unsigned int qp) const override;

protected:
  const bool _use_old_avg_prop;
  const MaterialProperty<Real> & _avg_mises_stress;
  const MaterialProperty<Real> & _avg_hyd_stress;
  const MaterialProperty<Real> & _avg_eq_strain;
  const MaterialProperty<Real> & _avg_eq_strain_rate;
  const MaterialProperty<std::vector<Real>> & _a_old;
  const MaterialProperty<std::vector<Real>> & _N_old;
  const MaterialProperty<std::vector<Real>> & _T_old;

  const Real _pi; /*convenint redeifinition of pi*/
  /// model paramters
  const Real _E_interface;   /*interface Young modulus*/
  const Real _E_penalty;     /*incase of copenetration*/
  const Real _beta_exponent; /*Traction exponent*/
  const Real _D_gb;          /*gran boundary diffusion coefficient*/
  const Real _n_exponent;    /*power law creep exponent*/
  const Real _alpha_n;       /*3/(2*_n_exponent))*/
  const Real _psi_angle;     /*equilibrium cavity tip half-angle [degree]*/
  const Real _h;             /*function of psi angle*/
  const Real _sigma_0;       /*traction normalization parameter*/
  const Real _S_thr;         /*threshold value of S for nucleation on set*/
  const Real _FN;            /*nucleation rate constant*/
  const Real _NI;            /*minimum cavity half spacing*/
  const Real _Nmax;          /*maximum cavity number density*/
  const Real _a0;            /*intial cavity radius*/
  const Real _b0;            /*minimum cavity half spacing*/
  const Real _b_sat;         /*minimum cavity half spacing*/

  /*var relate to non linear system for opening mode*/
  const std::vector<std::string> _var_names_newton_opening;
  const unsigned int _n_vars_newton_opening;
  const std::vector<std::string> _other_var_names_newton_opening;
  const unsigned int _n_other_var_names_newton_opening;
  const unsigned int _n_var_names_newton_opening_total;
  std::vector<std::string> _all_var_names_newton_opening;
  const bool _automatic_scale_opening_variable;
  const std::vector<Real> _standard_scale_factor;
  var4NL _opening_mode_vars;

  std::map<std::string, unsigned int> _map_var_name_idx_newton_opening;

  // DenseVector<Real> _newton_var_opening;
  // DenseVector<Real> _R_newton;
  // DenseMatrix<Real> _J_newton;
  // const bool _newton_scale_variables_opening;
  // const std::vector<Real> _standard_scale_factor_opening;
  //
  // std::map<std::string, unsigned int> _map_var_newton_idx;

  // methods for variable management during non linear calculation

  // newton general variables
  // std::vector<Real> _current_scale_factor;

private:
  // unsigned int getVarIndex(const std::string & /*var_name*/) const;
  // void zeroVarValues();
  // void assignVarValue(const std::string & /*var_name*/, const Real & /*value*/);
  // Real getVarValue(const std::string & /*var_name*/) const;
  //
  // void zeroDerivatives();
  // void assignDerivative(const std::string & /*var_name*/,
  //                       const std::string & /*dvar_name*/,
  //                       const Real & /*value*/);
  // Real getDerivative(const std::string & /*var_name*/, const std::string & /*dvar_name*/)
  // const;
  //
  // void
  // unrollVars(const DenseVector<Real> & current_x, Real & a, Real & b, Real & N, Real & TN)
  // const
  // {
  //   a = current_x(0);
  //   N = current_x(1);
  //   TN = current_x(2);
  //   b = fun_b_init(N);
  // }

  Real fun_b_init(const Real & N) const { return 1. / std::sqrt(N * _pi); }

  void fun_b_N();

  void computeNewtonResidualAndJacobian();

  void NewtonStep();

  void Newton();

  void computeCompliance();
  void Tn_fun();
  void TnElastic_Res();

  void N_Fun();
  void N_Res();

  void a_Fun();
  void a_Res();

  void dVi_Fun();
  void dV_1_Fun(Real & /*dv1*/, DenseVector<Real> & /*ddv1_dx*/, const unsigned int /*i_vol*/);
  void dV_2_Fun(Real & /*dv2*/, DenseVector<Real> & /*ddv2_dx*/, const unsigned int /*i_vol*/);
  void dVL2_Fun(Real & /*dvl2*/, DenseVector<Real> & /*ddvl2_dx*/);
  void dVH2_Fun(Real & /*dvh2*/, DenseVector<Real> & /*ddvh2_dx*/);
  void qi_Fun(Real & /*q*/, DenseVector<Real> & /*dq_dx*/, const unsigned int /*i_vol*/);
  void fi_Fun(Real & /*f*/, DenseVector<Real> & /*df_dx*/, const unsigned int /*i_vol*/);
  void fL_Fun(Real & /*fL*/, DenseVector<Real> & /*dfL_dx*/);
  void fH_Fun(Real & /*fH*/, DenseVector<Real> & /*dfH_dx*/);
  void f_a_over_b_square_Fun(Real & /*f_ab*/, DenseVector<Real> & /*df_ab_dx*/);
  void f_abL_Fun(Real & /*f_abL*/, DenseVector<Real> & /*df_abL_dx*/);
  int m_Fun();
  Real g_Fun();
  Real beta_n_m_Fun();

  Real ab_ratio(const bool & old = false);

  /// gives the total ec between the beginnin of Moose time step and current
  Real getDeltaEc();

  Real getEcRate();

  Real getAvgSVM();

  Real getAvgSH();

  /// gives the total duN between the beginnin of Moose time step and current
  Real getDeltauN();

  /// gives the total uN between the beginnin of Moose time step and current
  Real getuN();

  /// gives the  delta ec between the the previous newton solution and and current
  Real getDeltaEcSubStep();

  /// gives the  delta delta uN between the the previous newton solution and and current
  Real getDeltauNSubStep();

  /// gives the  delta ec between the the previous newton solution and and current
  Real getuNSubStep();

  void setZeroVectorDerivative(DenseVector<Real> & d_vector,
                               const std::string & varname,
                               const Real & value)
  {
    d_vector(_opening_mode_vars.getVarIndex(varname)) = value;
  };
};

#endif // CZMLAWSHAMNEEDLEMAN_H
