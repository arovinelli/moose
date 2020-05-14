//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceMaterial.h"

class CZMMaterialLD;
template <>
InputParameters validParams<CZMMaterialLD>();
/**
 * Implementation of the non-stateful exponential traction separation law
 * proposed by Salehani, Mohsen Khajeh and Irani, Nilgoon 2018
 **/
class CZMMaterialLD : public InterfaceMaterial
{
public:
  CZMMaterialLD(const InputParameters & parameters);

protected:
  virtual void initialSetup() override;

  // this method is responsible for computing the traction and traction rate in
  // the interface coordinate system
  virtual void computeTractionAndDerivatives();
  // virtual void computeTractionDerivatives();
  // virtual RankTwoTensor computeTractionDerivatives() override{};

  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// normal to the interface
  const MooseArray<Point> & _normals;

  /// number of displacement components
  const unsigned int _ndisp;

  /// the coupled displacement and neighbor displacement values
  ///@{
  std::vector<const VariableValue *> _disp;
  std::vector<const VariableValue *> _disp_neighbor;
  std::vector<MooseVariable *> _disp_vars;
  std::vector<const VariableValue *> _disp_old;
  std::vector<const VariableValue *> _disp_neighbor_old;
  ///@}

  /// the coupled test and trial fuinction
  ///@{
  std::vector<const VariablePhiValue *> _phi;
  std::vector<const VariablePhiValue *> _phi_neighbor;
  ///@}

  /// the coupled test and trial fuinction
  ///@{
  std::vector<const VariablePhiGradient *> _grad_phi;
  std::vector<const VariablePhiGradient *> _grad_phi_neighbor;
  ///@}

  std::vector<RankTwoTensor> _djump_global_duk;
  std::vector<RankTwoTensor> _djump_global_duk_neighbor;
  std::vector<RankThreeTensor> _dFmn_duk;
  std::vector<RankThreeTensor> _dFmn_neighbor_duk;

  /// the displacement jump in global and local coordiante
  ///@{
  MaterialProperty<RealVectorValue> & _displacement_jump_global;
  MaterialProperty<RealVectorValue> & _displacement_jump_global_old;
  MaterialProperty<RealVectorValue> & _displacement_jump;
  const MaterialProperty<RealVectorValue> & _displacement_jump_old;
  MaterialProperty<RealVectorValue> & _displacement_jump_global_inc;
  MaterialProperty<RealVectorValue> & _displacement_jump_inc;
  ///@}

  /// the value of the traction in global and local coordinates
  ///@{
  MaterialProperty<RealVectorValue> & _traction_global;
  MaterialProperty<RealVectorValue> & _traction_global_inc;
  MaterialProperty<RealVectorValue> & _traction;
  MaterialProperty<RealVectorValue> & _traction_inc;
  const MaterialProperty<RealVectorValue> & _traction_global_old;
  const MaterialProperty<RealVectorValue> & _traction_old;
  ///@}

  /// the traction's derivatives wrt the displacement jump in global and local
  /// coordinates
  ///@{
  MaterialProperty<std::vector<RankTwoTensor>> & _traction_derivatives_global;
  MaterialProperty<std::vector<RankTwoTensor>> & _traction_derivatives_global_neighbor;
  MaterialProperty<RankTwoTensor> & _traction_derivatives;
  ///@}

  // MaterialProperty<RealVectorValue> &_traction_global_deformed;
  // const MaterialProperty<RealVectorValue> &_traction_global_deformed_old;

  std::vector<const VariableGradient *> _grad_disp;
  std::vector<const VariableGradient *> _grad_disp_neighbor;

  // std::vector<const VariableGradient *> _grad_disp_inc;
  // std::vector<const VariableGradient *> _grad_disp_neighbor_inc;

  MaterialProperty<RankTwoTensor> & _F_avg;
  MaterialProperty<RankTwoTensor> & _DF_avg;
  const MaterialProperty<RankTwoTensor> & _F_avg_old;
  MaterialProperty<RankTwoTensor> & _R_avg;
  const MaterialProperty<RankTwoTensor> & _R_avg_old;
  MaterialProperty<RankTwoTensor> & _U_avg;
  MaterialProperty<RankTwoTensor> & _DR_avg;
  MaterialProperty<RankTwoTensor> & _DL_avg;
  const MaterialProperty<RankTwoTensor> & _DL_avg_old;
  MaterialProperty<RealVectorValue> & _n_avg;
  MaterialProperty<Real> & _dsdotdS_avg;
  MaterialProperty<Real> & _dadA_avg;

  // MaterialProperty<RankTwoTensor> &_linear_rot;
  // const MaterialProperty<RankTwoTensor> &_linear_rot_old;
  // MaterialProperty<RankTwoTensor> &_linear_rot_neighbor;
  // const MaterialProperty<RankTwoTensor> &_linear_rot_neighbor_old;
  // MaterialProperty<RankTwoTensor> &_linear_rot_avg;
  // const MaterialProperty<RankTwoTensor> &_linear_rot_avg_old;

  RankTwoTensor _K = RankTwoTensor::Identity();
  const bool _ld;
  const bool _use_area_change;

  // void update_L_df_vorticity(const RankTwoTensor &def_grad,
  //                            const RankTwoTensor &def_grad_old,
  //                            RankTwoTensor &L, RankTwoTensor &df,
  //                            RankTwoTensor &vorticity_inc);

  void computeFandL();

  void PolarDecomposition(const RankTwoTensor & F, RankTwoTensor & R, RankTwoTensor & U) const;

  // void computeRdot(const RankTwoTensor &F, const RankTwoTensor &L,
  //                  const RankTwoTensor &R, const RankTwoTensor &U,
  //                  RankTwoTensor &Rdot);

  void computeSdot(const RankTwoTensor & L, const RealVectorValue & n, Real & sdot);

  RankFourTensor computedRdF(const RankTwoTensor & R, const RankTwoTensor & U);

  RankTwoTensor RijkVk(const RankThreeTensor & R3, const RealVectorValue & V)
  {
    RankTwoTensor res = RankTwoTensor::initNone;
    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int j = 0; j < 3; j++)
        for (unsigned int k = 0; k < 3; k++)
          res(i, j) += R3(i, j, k) * V(k);
    return res;
  }

  RankTwoTensor RijklRkl(const RankFourTensor & R4, const RankTwoTensor & R2)
  {
    RankTwoTensor res = RankTwoTensor::initNone;
    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int j = 0; j < 3; j++)
        for (unsigned int k = 0; k < 3; k++)
          for (unsigned int l = 0; l < 3; l++)
            res(i, j) += R4(i, j, k, l) * R2(k, l);
    return res;
  }

  RankThreeTensor RjiklRklm(const RankFourTensor & R4, const RankThreeTensor & R3)
  {
    RankThreeTensor res;
    res.zero();
    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int j = 0; j < 3; j++)
        for (unsigned int k = 0; k < 3; k++)
          for (unsigned int l = 0; l < 3; l++)
            for (unsigned int m = 0; m < 3; m++)
              res(j, i, m) += R4(j, i, k, l) * R3(k, l, m);
    return res;
  }
};
