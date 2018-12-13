//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef CZMINTERFACEKERNELVISCOUS_H
#define CZMINTERFACEKERNELVISCOUS_H

#include "InterfaceTimeKernel.h"

/// Forward Declarations
class CZMInterfaceKernelViscous;

template <>
InputParameters validParams<CZMInterfaceKernelViscous>();

/// DG kernel implementing CZM for a 3D traction sepration law based on
/// the displacement jump. This kernel operates only on a single displacement
/// compenent. One kernel is needed for each dispalcement jump component
class CZMInterfaceKernelViscous : public InterfaceTimeKernel
{
public:
  CZMInterfaceKernelViscous(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type);
  virtual Real computeQpJacobian(Moose::DGJacobianType type);
  // virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar);

  /// the displacement component this kernel is operating on (0=x, 1=y, 2 =z)
  const unsigned int _disp_index;

  /// coupled displacement componenets values
  // const VariableValue & _disp_0_dot;
  // const VariableValue & _disp_0_dot_neighbor;
  const VariableValue & _disp_1_dot;
  const VariableValue & _disp_1_dot_neighbor;
  const VariableValue & _disp_2_dot;
  const VariableValue & _disp_2_dot_neighbor;

  /// coupled displacement componenets variables ID
  // unsigned int _disp_0_var;
  // unsigned int _disp_0_neighbor_var;
  unsigned int _disp_1_var;
  unsigned int _disp_1_neighbor_var;
  unsigned int _disp_2_var;
  unsigned int _disp_2_neighbor_var;

  /// variables containg the names of the material properties representing the
  /// reidual's and jacobian's coefficients. Derivates are assumed to be taken
  /// wrt to the displacement jump.
  const std::string _residual;
  const std::string _jacobian;

  const Real _viscosity_coefficient;
  // values of the residual's and jacobian's cofficients
  // const MaterialProperty<RealVectorValue> & _ResidualMP;
  // const MaterialProperty<std::vector<std::vector<Real>>> & _JacobianMP;
};

#endif // CZMINTERFACEKERNELVISCOUS_H
