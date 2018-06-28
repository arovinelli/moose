//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DisplacementJumpBasedCohesiveInterfaceMaterial.h"
#include "Assembly.h"
#include "MooseMesh.h"
// #include "TractionSeparationUOBase.h"
// #include "DisplacementJumpCohesiveInterface.h"
#include "RotationMatrix.h"

registerMooseObject("TensorMechanicsApp", DisplacementJumpBasedCohesiveInterfaceMaterial);

template <>
InputParameters
validParams<DisplacementJumpBasedCohesiveInterfaceMaterial>()
{
  InputParameters params = validParams<Material>();

  // params.addRequiredParam<UserObjectName>(
  //     "uo_TractionSeparationLaw",
  //     "the name of the user object including the traction separation law");
  params.addRequiredParam<UserObjectName>(
      "uo_CohesiveInterface", "the name of the user object implemnting the cohesive interface");
  params.addClassDescription("this material class is used when defining a "
                             "cohesive zone model to store stafeul properties");
  return params;
}

DisplacementJumpBasedCohesiveInterfaceMaterial::DisplacementJumpBasedCohesiveInterfaceMaterial(
    const InputParameters & parameters)
  : Material(parameters),

    _normals(_assembly.normals()),

    _uo_CohesiveInterface(getUserObject<DisplacementJumpCohesiveInterface>("uo_CohesiveInterface")),

    // _Jump(declareProperty<RealVectorValue>("Jump")),
    // _JumpLocal(declareProperty<RealVectorValue>("JumpLocal")),
    _Traction(declareProperty<RealVectorValue>("Traction")),
    _Traction_old(getMaterialPropertyOld<RealVectorValue>("Traction")),
    // _TractionLocal(&declareProperty<RealVectorValue>("TractionLocal")),
    _TractionSpatialDerivative(declareProperty<RankTwoTensor>("TractionSpatialDerivative"))
// _TractionSpatialDerivativeLocal(
//     &declareProperty<RankTwoTensor>("TractionSpatialDerivativeLocal")),

// _RotationGlobal2Local(RealTensorValue()),
// _RotationLocal2Global(RealTensorValue())

{
  // // assign user object
  // _uo_tractionSeparation = &getUserObjectByName<TractionSeparationUOBase>(
  //     parameters.get<UserObjectName>("uo_TractionSeparationLaw"));
  //
  // // get stateful material property number and names
  // _uo_tractionSeparation->statefulMaterialPropertyNames(_materialPropertyNames);
  // _num_stateful_material_properties = _materialPropertyNames.size();
  //
  // if (_num_stateful_material_properties > 0)
  // {
  //   // initialize the stateful material property values container
  //   _materialPropertyValues.resize(_num_stateful_material_properties);
  //   _materialPropertyValues_old.resize(_num_stateful_material_properties);
  //
  //   // declare properties
  //   for (unsigned int i = 0; i < _num_stateful_material_properties; ++i)
  //   {
  //     _materialPropertyValues[i] =
  //     &declareProperty<std::vector<Real>>(_materialPropertyNames[i]);
  //
  //     _materialPropertyValues_old[i] =
  //         &getMaterialPropertyOld<std::vector<Real>>(_materialPropertyNames[i]);
  //   }
  // }
}

void
DisplacementJumpBasedCohesiveInterfaceMaterial::computeQpProperties()
{

  RealVectorValue traction_temp;
  RankTwoTensor traction_spatial_derivative_temp;

  std::cout << "computeQpProperties" << std::endl;

  _uo_CohesiveInterface.computeResidaulAndJacobianCoefficients(
      _current_elem->id(), _current_side, _qp, traction_temp, traction_spatial_derivative_temp);

  std::cout << "save values" << std::endl;

  _Traction[_qp] = traction_temp;
  _TractionSpatialDerivative[_qp] = traction_spatial_derivative_temp;

  // // transform from global to loval cooridnates
  // moveToLocalFrame();
  //
  // // compute tractions and traction derivatives in localc coordinates
  // _uo_tractionSeparation->computeTractionLocal(_qp, (*_TractionLocal)[_qp]);
  //
  // _uo_tractionSeparation->computeTractionSpatialDerivativeLocal(
  //     _qp, (*_TractionSpatialDerivativeLocal)[_qp]);

  // // move results back to global coordinates
  // moveBackToGlobalFrame();
  std::cout << "finished computeQpProperties" << std::endl;
}

void
DisplacementJumpBasedCohesiveInterfaceMaterial::initQpStatefulProperties()
{

  for (unsigned int i = 0; i < 3; i++)
    _Traction[_qp](i) = -1;
}

// void
// DisplacementJumpBasedCohesiveInterfaceMaterial::moveToLocalFrame()
// {
//   // this is a rotation matrix that will rotate _n to the "x" axis such that the
//   // the first compoenent in the local frame represent the opening displacment
//   _RotationGlobal2Local = RotationMatrix::rotVec1ToVec2(_normals[_qp], RealVectorValue(1, 0, 0));
//
//   // compute the jump in the lcoal coordinate system
//   for (unsigned int i = 0; i < 3; i++)
//   {
//     _JumpLocal[_qp](i) = 0; // just to be sure that it reset at each iteration
//     for (unsigned int j = 0; j < 3; j++)
//     {
//       _JumpLocal[_qp](i) += _RotationGlobal2Local(i, j) * _Jump[_qp](j);
//     }
//   }
// }

// void
// DisplacementJumpBasedCohesiveInterfaceMaterial::moveBackToGlobalFrame()
// {
//
//   _RotationLocal2Global = _RotationGlobal2Local.transpose();
//   // rotate traction in the global frame
//   for (unsigned int i = 0; i < 3; i++)
//   {
//     _Traction[_qp](i) = 0;
//     for (unsigned int j = 0; j < 3; j++)
//     {
//       _Traction[_qp](i) += _RotationLocal2Global(i, j) * (*_TractionLocal)[_qp](j);
//     }
//   }
//
//   // rotate traction derivatives in the global frame
//   _TractionSpatialDerivative[_qp] = (*_TractionSpatialDerivativeLocal)[_qp];
//   _TractionSpatialDerivative[_qp].rotate(_RotationLocal2Global);
// }
