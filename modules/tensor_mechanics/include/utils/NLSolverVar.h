//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef NLSOLVERVAR_H
#define NLSOLVERVAR_H

#include "Moose.h"
#include "InputParameters.h"

// Any requisite includes here
#include "libmesh/dense_vector.h"
#include "libmesh/dense_matrix.h"

// class var4NL;
// template <>
// InputParameters validParams<var4NL>();

class var4NL
{
public:
  var4NL(std::vector<std::string> nls_var_names,
         std::vector<std::string> other_var_names,
         Real dt_newton,
         unsigned int qp);
  // var4NL();

  void setValue(const std::string & /*var_name*/, const Real & /*value*/);

  void setComputedValue(const std::string & /*var_name*/, const Real & /*value*/);

  void setResidual(const std::string & /*var_name*/, const Real & /*value*/);

  Real getValue(const std::string & /*var_name*/) const;

  Real getInitialValue(const std::string & /*var_name*/) const;

  Real getComputedValue(const std::string & /*var_name*/) const;

  Real getResidual(const std::string & /*var_name*/) const;

  Real getValueOld(const std::string & /*var_name*/) const;

  Real getValueScaled(const std::string & /*var_name*/);

  Real getValueScaledOld(const std::string & /*var_name*/);

  void setDerivative(const std::string & /*var_name*/,
                     const std::string & /*dvar_name*/,
                     const Real & /*value*/);

  Real getDerivative(const std::string & /*var_name*/, const std::string & /*dvar_name*/) const;

  DenseVector<Real> getDerivativeVector(const std::string & /*var_name*/);

  void setNLSScaleFactors();

  Real getScaleFactor(const std::string & /*var_name*/) const;

  DenseMatrix<Real> getNLSJacobian();
  void setNLSJacobian(const std::vector<std::string> & /*var_names*/,
                      const DenseMatrix<Real> & /*J*/);

  DenseMatrix<Real> getNLSJacobianScaled();

  DenseVector<Real> getNLSValue();

  DenseVector<Real> getNLSResidual();

  DenseVector<Real> getNLSValueScaled();

  DenseVector<Real> getNLSValueOld();

  DenseVector<Real> getNLSValueScaledOld();

  void setNLSValue(const std::vector<std::string> & /*var_names*/,
                   const DenseVector<Real> & /*var_values*/);

  void setNLSResidual(const std::vector<std::string> & /*var_names*/,
                      const DenseVector<Real> & /*var_residual*/);

  void initNLS(const std::vector<std::string> & /*var_names*/,
               const DenseVector<Real> & /*intial_values*/);

  void resetNLS();
  void advanceNLS() { _values_old = _values; };

  DenseVector<Real> getUnitDerivativeVector(const std::string & /*var_name*/);
  DenseVector<Real> getZeroDerivativeVector();

  Real getDt() { return _dt_newton; };
  void setDt(Real dt) { _dt_newton = dt; };
  void set_qp(unsigned int qp) { _qp = qp; };
  unsigned int get_qp() { return _qp; };
  unsigned int getVarIndex(const std::string & var_name) const;

protected:
  // const InputParameters & _params;
  void InitVarIdxMap();
  void zeroValues();
  void zeroDerivatives();
  void zeroResidual();

  Real _dt_newton;
  unsigned int _qp;
  const unsigned int _n_nls_vars;
  unsigned int _n_other_vars;
  unsigned int _n_vars_total;
  DenseVector<Real> _intial_values;
  const std::vector<std::string> _nls_var_names;
  const std::vector<std::string> _other_var_names;
  DenseVector<Real> _values;
  DenseVector<Real> _values_old;
  DenseVector<Real> _computed_values;
  DenseMatrix<Real> _derivatives;
  std::vector<Real> _current_scale_factor;
  DenseVector<Real> _var_NLS;
  DenseVector<Real> _R_NLS;
  DenseMatrix<Real> _Jac_NLS;
  std::map<std::string, unsigned int> _map_var_name_idx;
  std::vector<std::string> _all_var_names;
};

#endif // NLSOLVERVAR_H
