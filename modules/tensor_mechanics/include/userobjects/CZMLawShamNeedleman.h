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

  const Real _pi; /*convenint redeifinition of pi*/
  /// model paramters
  const Real _E_interface;   /*interface Young modulus*/
  const Real _beta_exponent; /*Traction exponent*/
  const Real _D_gb;          /*gran boundary diffusion coefficient*/
  const Real _n_exp;         /*power law creep exponent*/
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

  const unsigned int _n_vars_newton;
  DenseMatrix<Real> _J_newton;
  DenseVector<Real> _R_newton;
  DenseVector<Real> _X_newton;

private:
  Real fun_b(const Real N) { return (1. / std::sqrt(N * _pi)); }
  void NewtonStep(const std::vector<Real> /*x_old*/,
                  std::vector<Real> /*x_new*/,
                  std::vector<Real> /*residaul*/,
                  const std::vector<std::vector<Real>> /*J*/);
};

#endif // CZMLAWSHAMNEEDLEMAN_H
