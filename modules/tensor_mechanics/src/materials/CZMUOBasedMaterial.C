//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CZMUOBasedMaterial.h"
#include "RotationMatrix.h"

registerMooseObject("TensorMechanicsApp", CZMUOBasedMaterial);

template <>
InputParameters
validParams<CZMUOBasedMaterial>()
{
  InputParameters params = validParams<Material>();

  params.addRequiredParam<UserObjectName>(
      "traction_separation_UO",
      "the name of the user object including the traction separation law");
  // params.addParam<UserObjectName>(
  //     "unload_traction_separation_UO",
  //     "the name of the user object including the traction separation law");
  // params.addParam<UserObjectName>(
  //     "coopenetration_penalty_UO",
  //     "the name of the user object including the traction separation law");
  params.addRequiredParam<UserObjectName>(
      "displacement_jump_UO",
      "the name of the UO collecting bulk material property across an interface");
  params.addParam<Real>("coopenetration_penalty", 1, "copenation penalty factor");
  params.addClassDescription("this material class is used when defining a "
                             "cohesive zone model");
  params.addParam<bool>("compute_shear_traction", false, "flag to compute the shear traction");
  params.addParam<std::vector<std::string>>(
      "bulk_avg_mp_names", std::vector<std::string>(), "names of average mps");
  params.addParam<std::vector<UserObjectName>>(
      "bulk_avg_mp_uo_names", std::vector<UserObjectName>(0), "names of averaging uo");
  params.addParam<bool>(
      "need_avg_scalar_vars_derivatives", false, "flag to compute the shear traction");
  return params;
}

CZMUOBasedMaterial::CZMUOBasedMaterial(const InputParameters & parameters)
  : Material(parameters),
    _compute_shear_traction(getParam<bool>("compute_shear_traction")),
    _displacement_jump_UO(getUserObject<DispJumpUO_QP>("displacement_jump_UO")),
    _bulk_avg_mp_names(getParam<std::vector<std::string>>("bulk_avg_mp_names")),
    _n_bulk_avg_mp(getParam<std::vector<std::string>>("bulk_avg_mp_names").size()),
    _bulk_avg_mp_uo_names(getParam<std::vector<UserObjectName>>("bulk_avg_mp_uo_names")),
    _n_bulk_avg_mp_uo_names(getParam<std::vector<UserObjectName>>("bulk_avg_mp_uo_names").size()),
    _need_avg_scalar_vars_derivatives(getParam<bool>("need_avg_scalar_vars_derivatives")),
    _traction_derivatives_other_avg_vars_id(
        declareProperty<std::map<unsigned int /*var_id*/, unsigned int /*mp_index*/>>(
            "traction_derivatives_other_avg_vars_id")),
    _traction_separation_UO(getUserObject<CZMTractionSeparationUOBase>("traction_separation_UO")),
    // _unload_traction_separation_UO(
    //     parameters.isParamSetByUser("unload_traction_separation_UO")
    //         ? getUserObject<CZMTractionSeparationUOBase>("unload_traction_separation_UO")
    //         : _traction_separation_UO),
    // _coopenetration_penalty_UO(
    //     parameters.isParamSetByUser("coopenetration_penalty_UO")
    //         ? getUserObject<CZMTractionSeparationUOBase>("coopenetration_penalty_UO")
    //         : _traction_separation_UO),
    // _coopenetration_penalty(getParam<Real>("coopenetration_penalty")),
    // _selected_CZM_UO(&_unload_traction_separation_UO),
    _displacement_jump(declareProperty<RealVectorValue>("displacement_jump")),
    _displacement_jump_dot(declareProperty<RealVectorValue>("displacement_jump_dot")),
    _displacement_jump_local(declareProperty<RealVectorValue>("displacement_jump_local")),
    _displacement_jump_dot_local(declareProperty<RealVectorValue>("displacement_jump_dot_local")),
    _displacement_jump_local_old(
        getMaterialPropertyOld<RealVectorValue>("displacement_jump_local")),
    _traction(declareProperty<RealVectorValue>("traction")),
    _traction_local(declareProperty<RealVectorValue>("traction_local")),
    _traction_spatial_derivatives(declareProperty<RankTwoTensor>("traction_spatial_derivatives")),
    _traction_spatial_derivatives_local(
        declareProperty<RankTwoTensor>("traction_spatial_derivatives_local")),
    _czm_residual(declareProperty<RealVectorValue>("czm_residual")),
    _czm_jacobian(declareProperty<std::vector<std::vector<Real>>>("czm_jacobian")),
    _normals_average(declareProperty<RealVectorValue>("normals_average")),
    _shear_traction(declareProperty<Real>("shear_traction")),
    // _uo_id(0),
    _n_uo_czm_properties(_traction_separation_UO.getNumberStatefulMaterialProperties()),
    _n_non_stateful_uo_czm_properties(
        _traction_separation_UO.getNumberNonStatefulMaterialProperties())

{
  if (_n_uo_czm_properties > 0)
  {
    // initialize the userobject material property container
    _uo_czm_properties.resize(_n_uo_czm_properties);
    _uo_czm_properties_old.resize(_n_uo_czm_properties);
    for (unsigned int mp_index = 0; mp_index < _n_uo_czm_properties; mp_index++)
    {
      // declare a material property
      _uo_czm_properties[mp_index] = &declareProperty<std::vector<Real>>(
          _traction_separation_UO.getStatefulMaterialPropertyName(mp_index));
      // for a stateful material property get the old value
      _uo_czm_properties_old[mp_index] = &getMaterialPropertyOld<std::vector<Real>>(
          _traction_separation_UO.getStatefulMaterialPropertyName(mp_index));
    }
  }

  if (_n_non_stateful_uo_czm_properties > 0)
  {
    // initialize the userobject material property container
    _uo_non_stateful_czm_properties.resize(_n_non_stateful_uo_czm_properties);
    for (unsigned int mp_index = 0; mp_index < _n_non_stateful_uo_czm_properties; mp_index++)
    {
      // declare a material property
      _uo_non_stateful_czm_properties[mp_index] = &declareProperty<std::vector<Real>>(
          _traction_separation_UO.getNonStatefulMaterialPropertyName(mp_index));
    }
  }

  if (_n_bulk_avg_mp != _n_bulk_avg_mp_uo_names)
    mooseError("CZMUOBasedMaterial:: the number fo required averaged material properties does not "
               "match the number of provided averaging user objesct");

  // assgin user objects and declare associate material properties
  if (_n_bulk_avg_mp_uo_names > 0)
  {
    // assign user obejcts
    _bulk_avg_mp_uo.resize(_n_bulk_avg_mp_uo_names);
    _avg_interface_mp.resize(_n_bulk_avg_mp_uo_names);
    _traction_derivatives_other_avg_vars_local.resize(_n_bulk_avg_mp_uo_names);
    _traction_derivatives_other_avg_vars.resize(_n_bulk_avg_mp_uo_names);
    for (unsigned int i = 0; i < _n_bulk_avg_mp_uo_names; ++i)
    {
      _bulk_avg_mp_uo[i] =
          &getUserObjectByName<ScalarBulkMPAcrossInterface_QP>(_bulk_avg_mp_uo_names[i]);
      _avg_interface_mp[i] = &declareProperty<Real>(_bulk_avg_mp_names[i]);
      if (_need_avg_scalar_vars_derivatives)
      {
        _traction_derivatives_other_avg_vars_local[i] =
            &declareProperty<RealVectorValue>(_bulk_avg_mp_names[i] + "_der_local");
        _traction_derivatives_other_avg_vars[i] =
            &declareProperty<RealVectorValue>(_bulk_avg_mp_names[i] + "_der");
      }
    }
  }
}

void
CZMUOBasedMaterial::computeQpProperties()
{
  // update average interface properties and save varaibles id
  _traction_derivatives_other_avg_vars_id[_qp].clear();

  if (_n_bulk_avg_mp_uo_names > 0)
    for (unsigned int avg_mp_index = 0; avg_mp_index < _n_bulk_avg_mp_uo_names; avg_mp_index++)
    {
      (*_avg_interface_mp[avg_mp_index])[_qp] =
          _bulk_avg_mp_uo[avg_mp_index]->getValueAverage(_current_elem->id(), _current_side, _qp);
      _traction_derivatives_other_avg_vars_id[_qp]
                                             [_bulk_avg_mp_uo[avg_mp_index]->getScalarVarID()] =
                                                 avg_mp_index + 3;
    }
  // resize non stateful mp
  if (_n_non_stateful_uo_czm_properties > 0)
    for (unsigned int mp_index = 0; mp_index < _n_non_stateful_uo_czm_properties; mp_index++)
      (*_uo_non_stateful_czm_properties[mp_index])[_qp].resize(
          _traction_separation_UO.getNonStatefulMaterialPropertySize(mp_index));

  for (unsigned int i = 0; i < 3; i++)
  {
    // _normals_average[_qp](i) = (_normals_MP[_qp](i) + _normals_neighbor_MP[_qp](i)) / 2;
    _normals_average[_qp](i) = _normals[_qp](i);
  }
  RealTensorValue RotationGlobal2Local =
      RotationMatrix::rotVec1ToVec2(_normals_average[_qp], RealVectorValue(1, 0, 0));

  _displacement_jump[_qp] =
      _displacement_jump_UO.getDisplacementJump(_current_elem->id(), _current_side, _qp);

  _displacement_jump_local[_qp] = rotateVector(_displacement_jump[_qp], RotationGlobal2Local);

  _displacement_jump_dot[_qp] =
      _displacement_jump_UO.getDisplacementJumpVelocity(_current_elem->id(), _current_side, _qp);

  _displacement_jump_dot_local[_qp] =
      rotateVector(_displacement_jump_dot[_qp], RotationGlobal2Local);

  if (_n_uo_czm_properties > 0)
    for (unsigned int mp_index = 0; mp_index < _n_uo_czm_properties; mp_index++)
      (*_uo_czm_properties[mp_index])[_qp] =
          _traction_separation_UO.getNewStatefulMaterialProperty(_qp, mp_index);
  if (_n_non_stateful_uo_czm_properties > 0)
    for (unsigned int mp_index = 0; mp_index < _n_non_stateful_uo_czm_properties; mp_index++)
      (*_uo_non_stateful_czm_properties[mp_index])[_qp] =
          _traction_separation_UO.getNewNonStatefulMaterialProperty(_qp, mp_index);

  _traction_local[_qp] = _traction_separation_UO.computeTractionLocal(_qp);
  //
  _traction_spatial_derivatives_local[_qp] =
      _traction_separation_UO.computeTractionSpatialDerivativeLocal(_qp);

  _traction[_qp] = rotateVector(_traction_local[_qp], RotationGlobal2Local, /*inverse =*/true);
  _traction_spatial_derivatives[_qp] = rotateTensor2(
      _traction_spatial_derivatives_local[_qp], RotationGlobal2Local, /*inverse =*/true);

  _czm_jacobian[_qp].resize(3, std::vector<Real>(3, 0));

  if (_n_bulk_avg_mp_uo_names > 0 && _need_avg_scalar_vars_derivatives)
  {
    _czm_jacobian[_qp].resize(3 + _n_bulk_avg_mp_uo_names, std::vector<Real>(3, 0));
    for (unsigned int avg_mp_index = 0; avg_mp_index < _n_bulk_avg_mp_uo_names; avg_mp_index++)
    {
      (*_traction_derivatives_other_avg_vars_local[avg_mp_index])[_qp] =
          _traction_separation_UO.computeTractionOtherAveragedScalarVarDerivatives(_qp,
                                                                                   avg_mp_index);
      (*_traction_derivatives_other_avg_vars[avg_mp_index])[_qp] =
          rotateVector((*_traction_derivatives_other_avg_vars_local[avg_mp_index])[_qp],
                       RotationGlobal2Local,
                       /*inverse =*/true);
    }
  }

  _czm_residual[_qp] = _traction[_qp];
  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      _czm_jacobian[_qp][i][j] = _traction_spatial_derivatives[_qp](i, j);

  if (_n_bulk_avg_mp_uo_names > 0 && _need_avg_scalar_vars_derivatives)
    for (unsigned int i = 3; i < 3 + _n_bulk_avg_mp_uo_names; i++)
      for (unsigned int j = 0; j < 3; j++)
        _czm_jacobian[_qp][i][j] = (*_traction_derivatives_other_avg_vars[i - 3])[_qp](j);

  if (_compute_shear_traction)
    _shear_traction[_qp] =
        std::sqrt(std::pow(_traction_local[_qp](1), 2) + std::pow(_traction_local[_qp](2), 2));
}

void
CZMUOBasedMaterial::initQpStatefulProperties()
{
  RealTensorValue RotationGlobal2Local =
      RotationMatrix::rotVec1ToVec2(_normals[_qp], RealVectorValue(1, 0, 0));

  for (unsigned int i = 0; i < 3; i++)
  {
    _displacement_jump[_qp](i) = 0;
    _displacement_jump_local[_qp](i) = 0;
  }

  /// is there a special procedure for reload?
  if (_n_uo_czm_properties > 0)
  {
    for (unsigned int mp_index = 0; mp_index < _n_uo_czm_properties; mp_index++)
    {
      (*_uo_czm_properties[mp_index])[_qp].resize(
          _traction_separation_UO.getStatefulMaterialPropertySize(mp_index));

      (*_uo_czm_properties[mp_index])[_qp] =
          _traction_separation_UO.getStatefulMaterialPropertysIntialValues(mp_index);
    }
  }
}

RealVectorValue
CZMUOBasedMaterial::rotateVector(const RealVectorValue v,
                                 const RealTensorValue R,
                                 const bool inverse /*= false*/)
{
  RealTensorValue R_loc = R;
  if (inverse)
    R_loc = R_loc.transpose();

  RealVectorValue vrot;

  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = 0; j < 3; j++)
      vrot(i) += v(j) * R_loc(i, j);
  return vrot;
}

RankTwoTensor
CZMUOBasedMaterial::rotateTensor2(const RankTwoTensor T,
                                  const RealTensorValue R,
                                  const bool inverse /*= false*/)
{
  RealTensorValue R_loc = R;
  if (inverse)
    R_loc = R_loc.transpose();

  RankTwoTensor trot = T;
  trot.rotate(R_loc);
  return trot;
}
