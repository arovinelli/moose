//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef DISPJUMPACROSSINTERFACE_H
#define DISPJUMPACROSSINTERFACE_H

#include "InterfaceUserObject.h"

class DispJumpAcrossInterface;

template <>
InputParameters validParams<DispJumpAcrossInterface>();

class DispJumpAcrossInterface : public InterfaceUserObject
{
public:
  DispJumpAcrossInterface(const InputParameters & parameters);
  virtual ~DispJumpAcrossInterface();

  virtual RealVectorValue getDisplacementJump(const dof_id_type /*elem_id*/,
                                              const dof_id_type /*side_id*/,
                                              const unsigned int /*qp*/);

  virtual void initialize();
  virtual void execute();
  virtual void finalize();
  virtual void threadJoin(const UserObject & uo);

protected:
  std::map<std::pair<dof_id_type, unsigned int>, std::vector<std::vector<RealVectorValue>>>
      _map_values;

  const VariableValue & _disp_x;
  const VariableValue & _disp_x_neighbor;
  const VariableValue & _disp_y;
  const VariableValue & _disp_y_neighbor;
  const VariableValue & _disp_z;
  const VariableValue & _disp_z_neighbor;

  virtual RealVectorValue computeDisplacementJump(unsigned int /*qp*/);
};

#endif // DISPJUMPACROSSINTERFACE_H
