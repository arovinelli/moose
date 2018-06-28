//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef DISPLACEMENTJUMPCOHESIVEINTERFACE_H
#define DISPLACEMENTJUMPCOHESIVEINTERFACE_H

#include "InternalSideUserObject.h"

class DisplacementJumpCohesiveInterface;

template <>
InputParameters validParams<DisplacementJumpCohesiveInterface>();

class DisplacementJumpCohesiveInterface : public InternalSideUserObject
{
public:
  DisplacementJumpCohesiveInterface(const InputParameters & parameters);
  virtual ~DisplacementJumpCohesiveInterface();

  virtual void initialize();
  virtual void execute();
  virtual void finalize();
  virtual void threadJoin(const UserObject & uo);

  /// @{ Block all methods that are not used in explicitly called UOs
  // virtual void initialize() override;
  // virtual void execute() override;
  // virtual void finalize() override;
  // virtual void threadJoin(const UserObject &) override;

  void computeResidaulAndJacobianCoefficients(dof_id_type /*elem*/,
                                              unsigned int /*side*/,
                                              unsigned int /*qp*/,
                                              RealVectorValue & /*Traction*/,
                                              RankTwoTensor & /*TractionSpatialDerivative*/) const;

protected:
  std::map<std::pair<dof_id_type, unsigned int>, std::vector<std::vector<RealVectorValue>>>
      _map_values;

  const VariableValue & _disp_x;
  const VariableValue & _disp_x_neighbor;
  const VariableValue & _disp_y;
  const VariableValue & _disp_y_neighbor;
  const VariableValue & _disp_z;
  const VariableValue & _disp_z_neighbor;

  void moveToLocalFrame(RealVectorValue & /*Jump*/,
                        RealVectorValue & /*JumpLocal*/,
                        RealVectorValue & /*normals*/,
                        RealTensorValue & /*RotationGlobal2Local*/) const;

  void moveBackToGlobalFrame(RealVectorValue & /*TractionLocal*/,
                             RankTwoTensor & /*TractionSpatialDerivativeLocal*/,
                             RealVectorValue & /*Traction*/,
                             RankTwoTensor & /*TractionSpatialDerivative*/,
                             RealTensorValue & /*RotationGlobal2Local*/) const;

  virtual RealVectorValue computeDisplacementJump(unsigned int /*qp*/);

  virtual void computeTractionLocal(RealVectorValue & /*_JumpLocal*/,
                                    RealVectorValue & /*TractionLocal*/) const;

  virtual void
  computeTractionSpatialDerivativeLocal(RealVectorValue & /*_JumpLocal*/,
                                        RankTwoTensor & /*TractionDerivativeLocal*/) const;

private:
  // cohesive law paramters
  const std::vector<Real> _deltaU0;
  const std::vector<Real> _maxAllowableTraction;
};

#endif // DISPLACEMENTJUMPCOHESIVEINTERFACE_H
