//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef CZMUOBASEDMATERIAL_H
#define CZMUOBASEDMATERIAL_H

#include "Material.h"
#include "DispJumpUO_QP.h"
#include "ScalarBulkMPAcrossInterface_QP.h"
#include "CZMTractionSeparationUOBase.h"
class CZMUOBasedMaterial;
template <>
InputParameters validParams<CZMUOBasedMaterial>();
/**
 *
 */
class CZMUOBasedMaterial : public Material
{
public:
  CZMUOBasedMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;

  const bool _compute_shear_traction;

  /// User objects computing the displacement jump
  const DispJumpUO_QP & _displacement_jump_UO;

  /// name of require averaged bulk mp
  const std::vector<std::string> _bulk_avg_mp_names;
  /// number of require averaged bulk mp (inherited from _bulk_avg_mp_names length)
  const unsigned int _n_bulk_avg_mp;
  /// userobjects returning the require bilk mp averaged
  const std::vector<UserObjectName> _bulk_avg_mp_uo_names;
  const unsigned int _n_bulk_avg_mp_uo_names;
  std::vector<const ScalarBulkMPAcrossInterface_QP *> _bulk_avg_mp_uo;
  std::vector<MaterialProperty<Real> *> _avg_interface_mp;
  std::vector<MaterialProperty<RealVectorValue> *> _traction_derivatives_other_avg_vars;
  std::vector<MaterialProperty<RealVectorValue> *> _traction_derivatives_other_avg_vars_local;
  const bool _need_avg_scalar_vars_derivatives;
  MaterialProperty<std::map<unsigned int /*var_id*/, unsigned int /*mp_index*/>> &
      _traction_derivatives_other_avg_vars_id;

  /// User objectets defining the traction separation law
  /// non linear TS law
  const CZMTractionSeparationUOBase & _traction_separation_UO;
  // /// unloading behavior for stateful laws
  // const CZMTractionSeparationUOBase & _unload_traction_separation_UO;
  // /// penalty for copentration behavior
  // const CZMTractionSeparationUOBase & _coopenetration_penalty_UO;

  // const Real _coopenetration_penalty;

  /// penalty for copentration behavior
  // const CZMTractionSeparationUOBase * _selected_CZM_UO;

  /// the dispalcement jump in global coordiantes
  MaterialProperty<RealVectorValue> & _displacement_jump;

  /// the dispalcement jump in global coordiantes
  MaterialProperty<RealVectorValue> & _displacement_jump_dot;

  /// the disaplcement jump in natural element coordiantes
  MaterialProperty<RealVectorValue> & _displacement_jump_local;

  /// the disaplcement jump in natural element coordiantes
  MaterialProperty<RealVectorValue> & _displacement_jump_dot_local;

  /// the disaplcement jump in natural element coordiantes at the previous time step
  const MaterialProperty<RealVectorValue> & _displacement_jump_local_old;

  /// the value of the Traction in global coordiantes
  MaterialProperty<RealVectorValue> & _traction;

  /// the value of the Traction in natural element coordiantes
  MaterialProperty<RealVectorValue> & _traction_local;

  /// the value of the traction derivatives in global coordiantes
  MaterialProperty<RankTwoTensor> & _traction_spatial_derivatives;

  /// the value of the traction derivatives in natural element coordiantes
  MaterialProperty<RankTwoTensor> & _traction_spatial_derivatives_local;

  /// the material property in which the residual is stored
  MaterialProperty<RealVectorValue> & _czm_residual;

  /// the material property in which the jacobian is stored
  MaterialProperty<std::vector<std::vector<Real>>> & _czm_jacobian;

  // MaterialProperty<RealVectorValue> & _normals_MP;
  // MaterialProperty<RealVectorValue> & _normals_neighbor_MP;
  MaterialProperty<RealVectorValue> & _normals_average;

  MaterialProperty<Real> & _shear_traction;

  // unsigned int _uo_id;

  const unsigned int _n_uo_czm_properties;
  const unsigned int _n_non_stateful_uo_czm_properties;
  std::vector<MaterialProperty<std::vector<Real>> *> _uo_czm_properties;
  std::vector<const MaterialProperty<std::vector<Real>> *> _uo_czm_properties_old;
  std::vector<MaterialProperty<std::vector<Real>> *> _uo_non_stateful_czm_properties;

  // /// method selecting proper UO depending on laoding state
  // void selectCzmUO();

  /// Rotate a vector "T" via the rotation matrix "R".
  /// inverse rotation is achieved by setting "inverse" = true
  RealVectorValue rotateVector(const RealVectorValue /*V*/,
                               const RealTensorValue /*R*/,
                               const bool inverse = false);

  /// Rotate a rank2 tensor "T" via the rotation matrix "R".
  /// inverse rotation is achieved by setting "inverse" = true
  RankTwoTensor
  rotateTensor2(const RankTwoTensor /*T*/, const RealTensorValue /*R*/, const bool inverse = false);
};

#endif // CZMUOBASEDMATERIAL_H
