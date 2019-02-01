//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NLSolverVar.h"
#include "MooseError.h"

// template <>
// InputParameters
// validParams<var4NL>()
// {
//   InputParameters _params = emptyInputParameters();
//   _params.addRequiredParam<std::vector<std::string>>("nls_var_names",
//                                                      "name of the nonlinear variables");
//   _params.addRequiredParam<std::vector<Real>>("nls_var_values", "name of the nonlinear
//   variables"); _params.addRequiredParam<std::vector<std::string>>("other_var_names", "name of
//   other variables"); return _params;
// }
var4NL::var4NL(std::vector<std::string> nls_var_names,
               std::vector<std::string> other_var_names,
               Real dt_newton,
               unsigned int qp)
  : _dt_newton(dt_newton),
    _qp(qp),
    _n_nls_vars(nls_var_names.size()),
    _n_other_vars(other_var_names.size()),
    _n_vars_total(_n_nls_vars + _n_other_vars),
    _intial_values(_n_vars_total),
    _nls_var_names(nls_var_names),
    _other_var_names(other_var_names),
    _values(_n_vars_total),
    _values_old(_n_vars_total),
    _computed_values(_n_vars_total),
    _derivatives(_n_vars_total, _n_vars_total),
    _current_scale_factor(_n_nls_vars),
    _var_NLS(_n_nls_vars),
    _R_NLS(_n_nls_vars),
    _Jac_NLS(_n_nls_vars, _n_nls_vars),
    _all_var_names(_n_vars_total)
{
  InitVarIdxMap();
  // zero var and derivatives
  zeroValues();
  zeroDerivatives();
  zeroResidual();
}

void
var4NL::InitVarIdxMap()
{
  _map_var_name_idx.clear();
  // initialize varname_index map
  for (unsigned int i = 0; i < _n_nls_vars; i++)
  {
    _map_var_name_idx[_nls_var_names[i]] = i;
    _all_var_names[i] = _nls_var_names[i];
  }

  for (unsigned int i = _n_nls_vars; i < _n_vars_total; i++)
  {
    _map_var_name_idx[_other_var_names[i]] = i;
    _all_var_names[i] = _other_var_names[i];
  }
}

unsigned int
var4NL::getVarIndex(const std::string & var_name) const
{
  auto var_idx_p = _map_var_name_idx.find(var_name);
  if (var_idx_p != _map_var_name_idx.end())
    return var_idx_p->second;
  else
    mooseError("var4NL::getVarIndex cannot find variable named ", var_name);
}

void
var4NL::zeroValues()
{
  for (unsigned int i = 0; i < _n_nls_vars; i++)
  {
    _values(i) = 0.;
    _values_old(i) = 0.;
    _computed_values(i) = 0.;
  }
}

void
var4NL::zeroResidual()
{
  for (unsigned int i = 0; i < _n_nls_vars; i++)
    _R_NLS(i) = 0.;
}

void
var4NL::setValue(const std::string & var_name, const Real & value)
{
  unsigned int index_var = getVarIndex(var_name);
  _values(index_var) = value;
}

void
var4NL::setComputedValue(const std::string & var_name, const Real & value)
{
  unsigned int index_var = getVarIndex(var_name);
  _computed_values(index_var) = value;
}

void
var4NL::setResidual(const std::string & var_name, const Real & value)
{
  unsigned int index_var = getVarIndex(var_name);
  _R_NLS(index_var) = value;
}

Real
var4NL::getValue(const std::string & var_name) const
{
  unsigned int index_var = getVarIndex(var_name);
  return _values(index_var);
}

Real
var4NL::getInitialValue(const std::string & var_name) const
{
  unsigned int index_var = getVarIndex(var_name);
  return _intial_values(index_var);
}

Real
var4NL::getComputedValue(const std::string & var_name) const
{
  unsigned int index_var = getVarIndex(var_name);
  return _computed_values(index_var);
}

Real
var4NL::getResidual(const std::string & var_name) const
{
  unsigned int index_var = getVarIndex(var_name);
  return _R_NLS(index_var);
}

Real
var4NL::getValueOld(const std::string & var_name) const
{
  unsigned int index_var = getVarIndex(var_name);
  return _values_old(index_var);
}

Real
var4NL::getValueScaled(const std::string & var_name)
{
  return getValue(var_name) / getScaleFactor(var_name);
}

Real
var4NL::getValueScaledOld(const std::string & var_name)
{
  return getValueOld(var_name) / getScaleFactor(var_name);
}

void
var4NL::zeroDerivatives()
{
  for (unsigned int i = 0; i < _n_vars_total; i++)
    for (unsigned int j = 0; j < _n_vars_total; j++)
      _derivatives(i, j) = 0.;
}

void
var4NL::setDerivative(const std::string & var_name,
                      const std::string & dvar_name,
                      const Real & value)
{
  unsigned int index_var = getVarIndex(var_name);
  unsigned int index_dvar = getVarIndex(dvar_name);

  _derivatives(index_var, index_dvar) = value;
}

Real
var4NL::getDerivative(const std::string & var_name, const std::string & dvar_name) const
{
  return _derivatives(getVarIndex(var_name), getVarIndex(dvar_name));
}

DenseVector<Real>
var4NL::getDerivativeVector(const std::string & var_name)
{
  DenseVector<Real> d_vector(_n_vars_total);
  unsigned int var_idx = getVarIndex(var_name);
  for (unsigned int i = 0; i < _n_vars_total; i++)
    d_vector(i) = _derivatives(var_idx, i);
  return d_vector;
}

void
var4NL::setNLSScaleFactors()
{
  for (unsigned int i = 0; i < _n_nls_vars; i++)
    // if (std::abs(_values(i)) > 0)
    //   _current_scale_factor[i] = std::abs(_values(i));
    // else
    _current_scale_factor[i] = 1.;
}

Real
var4NL::getScaleFactor(const std::string & var_name) const
{
  unsigned int v_idx = getVarIndex(var_name);
  if (v_idx < _n_nls_vars)
    return _current_scale_factor[v_idx];
  else
    return 1;
}

DenseMatrix<Real>
var4NL::getNLSJacobian()
{
  for (unsigned int i = 0; i < _n_nls_vars; i++)
    for (unsigned int j = 0; j < _n_nls_vars; j++)
      _Jac_NLS(i, j) = getDerivative(_nls_var_names[i], _nls_var_names[j]);

  return _Jac_NLS;
}

void
var4NL::setNLSJacobian(const std::vector<std::string> & var_names, const DenseMatrix<Real> & J)
{
  for (unsigned int i = 0; i < _n_nls_vars; i++)
    for (unsigned int j = 0; j < _n_nls_vars; j++)
      setDerivative(var_names[i], var_names[j], J(i, j));
}

DenseMatrix<Real>
var4NL::getNLSJacobianScaled()
{
  _Jac_NLS = getNLSJacobian();
  for (unsigned int i = 0; i < _n_nls_vars; i++)
    for (unsigned int j = 0; j < _n_nls_vars; j++)
      _Jac_NLS(i, j) /= _current_scale_factor[j];

  return _Jac_NLS;
}

DenseVector<Real>
var4NL::getNLSValue()
{
  for (unsigned int i = 0; i < _n_nls_vars; i++)
    _var_NLS(i) = getValue(_nls_var_names[i]);
  return _var_NLS;
}

DenseVector<Real>
var4NL::getNLSResidual()
{
  return _R_NLS;
}

DenseVector<Real>
var4NL::getNLSValueScaled()
{
  _var_NLS = getNLSValue();
  for (unsigned int i = 0; i < _n_nls_vars; i++)
    _var_NLS(i) /= _current_scale_factor[i];
  return _var_NLS;
}

DenseVector<Real>
var4NL::getNLSValueOld()
{
  for (unsigned int i = 0; i < _n_nls_vars; i++)
    _var_NLS(i) = getValueOld(_nls_var_names[i]);
  return _var_NLS;
}

DenseVector<Real>
var4NL::getNLSValueScaledOld()
{
  _var_NLS = getNLSValueOld();
  for (unsigned int i = 0; i < _n_nls_vars; i++)
    _var_NLS(i) /= _current_scale_factor[i];
  return _var_NLS;
}

void
var4NL::setNLSValue(const std::vector<std::string> & var_names,
                    const DenseVector<Real> & var_values)
{
  for (unsigned int i = 0; i < _n_nls_vars; i++)
    setValue(var_names[i], var_values(i));
}

void
var4NL::setNLSResidual(const std::vector<std::string> & var_names,
                       const DenseVector<Real> & var_residual)
{
  for (unsigned int i = 0; i < _n_nls_vars; i++)
    setResidual(var_names[i], var_residual(i));
}

void
var4NL::initNLS(const std::vector<std::string> & var_names, const DenseVector<Real> & intial_values)
{
  for (unsigned int i = 0; i < _n_vars_total; i++)
  {
    unsigned int var_idx = getVarIndex(var_names[i]);
    _intial_values(var_idx) = intial_values(i);
    _values(var_idx) = intial_values(i);
    _values_old(var_idx) = intial_values(i);
    _R_NLS(var_idx) = 0;
  }
  setNLSScaleFactors();
}

void
var4NL::resetNLS()
{
  zeroValues();
  zeroDerivatives();
  zeroResidual();
  setNLSValue(_nls_var_names, _intial_values);
}

DenseVector<Real>
var4NL::getUnitDerivativeVector(const std::string & var_name)
{
  DenseVector<Real> unitPD = getZeroDerivativeVector();
  unsigned int var_idx = getVarIndex(var_name);
  unitPD(var_idx) = 1;
  return unitPD;
}

DenseVector<Real>
var4NL::getZeroDerivativeVector()
{
  DenseVector<Real> zeroPD(_n_vars_total);
  for (unsigned int i = 0; i < _n_vars_total; i++)
    zeroPD(i) = 0;
  return zeroPD;
}
