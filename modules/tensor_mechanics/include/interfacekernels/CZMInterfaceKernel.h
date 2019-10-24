//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceKernel.h"

/// Forward Declarations
class CZMInterfaceKernel;

template <>
InputParameters validParams<CZMInterfaceKernel>();

/// DG kernel implementing CZM for a 3D traction sepration law based on
/// the displacement jump. This kernel operates only on a single displacement
/// compenent. One kernel is needed for each dispalcement jump component
class CZMInterfaceKernel : public InterfaceKernel
{
public:
  CZMInterfaceKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type);
  virtual Real computeQpJacobian(Moose::DGJacobianType type);
  virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar);

  /// the displacement component this kernel is operating on (0=x, 1=y, 2 =z)
  const unsigned int _disp_index;

  /// number of displacement components
  const unsigned int _ndisp;

  /// coupled displacement componenets variables ID
  std::vector<unsigned int> _disp_var;
  std::vector<unsigned int> _disp_neighbor_var;

  // values of the residual's and jacobian's cofficients
  const MaterialProperty<RealVectorValue> & _traction;
  const MaterialProperty<RankTwoTensor> & _traction_derivative;
};
