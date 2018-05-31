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
#include "BoundaryRestrictableRequired.h"

class DisplacementJumpCohesiveInterface;

template <>
InputParameters validParams<DisplacementJumpCohesiveInterface>();

class DisplacementJumpCohesiveInterface : public InternalSideUserObject,
                                          public BoundaryRestrictableRequired
{
public:
  DisplacementJumpCohesiveInterface(const InputParameters & parameters);

  /// @{ Block all methods that are not used in explicitly called UOs
  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;
  virtual void threadJoin(const UserObject &) override final;

  std::map<std::pair<dof_id_type, unsigned int>, std::vector<RealVectorValue>> _map_values;

protected:
  const VariableValue & _disp_x;
  const VariableValue & _disp_x_neighbor;
  const VariableValue & _disp_y;
  const VariableValue & _disp_y_neighbor;
  const VariableValue & _disp_z;
  const VariableValue & _disp_z_neighbor;

private:
  virtual RealVectorValue computeDisplacementJump();
};

#endif // DISPLACEMENTJUMPCOHESIVEINTERFACE_H
