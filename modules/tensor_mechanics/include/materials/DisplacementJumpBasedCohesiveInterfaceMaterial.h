//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef DISPLACEMENTJUMPBASEDCOHESIVEINTERFACEMATERIAL_H
#define DISPLACEMENTJUMPBASEDCOHESIVEINTERFACEMATERIAL_H

#include "Material.h"
// #include "TractionSeparationUOBase.h"
#include "DisplacementJumpCohesiveInterface.h"

class DisplacementJumpBasedCohesiveInterfaceMaterial;

template <>
InputParameters validParams<DisplacementJumpBasedCohesiveInterfaceMaterial>();

/**
 *
 */
class DisplacementJumpBasedCohesiveInterfaceMaterial : public Material
{
public:
  DisplacementJumpBasedCohesiveInterfaceMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;

  const MooseArray<Point> & _normals;

  /// User objects defining the traction separation law
  // const TractionSeparationUOBase * _uo_tractionSeparation;

  /// User objects computing the displacement jump
  const DisplacementJumpCohesiveInterface & _uo_CohesiveInterface;

  // std::vector<std::string> _materialPropertyNames;
  // unsigned int _num_stateful_material_properties;
  //
  // std::vector<MaterialProperty<std::vector<Real>> *> _materialPropertyValues;
  // std::vector<const MaterialProperty<std::vector<Real>> *> _materialPropertyValues_old;

  /// the disaplcement jump in global coordiantes
  // MaterialProperty<RealVectorValue> & _Jump;
  //
  // /// the disaplcement jump in natural element coordiantes
  // MaterialProperty<RealVectorValue> & _JumpLocal;

  /// the value of the Traction in global coordiantes
  MaterialProperty<RealVectorValue> & _Traction;
  const MaterialProperty<RealVectorValue> & _Traction_old;

  /// the value of the Traction in natural element coordiantes
  // MaterialProperty<RealVectorValue> * _TractionLocal;

  /// the value of the Traction Derivatives in global coordiantes
  MaterialProperty<RankTwoTensor> & _TractionSpatialDerivative;

  /// the value of the Traction Derivatives in natural element coordiantes
  // MaterialProperty<RankTwoTensor> * _TractionSpatialDerivativeLocal;

  /// rotation matrix that takes _n to (0, 0, 1) i.e. _nLocal = R*_n
  /// global to local and viceversa
  // RealTensorValue _RotationGlobal2Local;
  // RealTensorValue _RotationLocal2Global;

  /// compute the displacmenet Jump in natural element coordinates (N,T,S)
  // virtual void moveToLocalFrame();

  /// Transform Travtion and Traction derivatives from the local back to global
  // virtual void moveBackToGlobalFrame();
};

#endif // DISPLACEMENTJUMPBASEDCOHESIVEINTERFACEMATERIAL_H
