//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "InterfaceValueUserObject.h"

class InterfaceQpValueUserObjectBase;

template <>
InputParameters validParams<InterfaceQpValueUserObjectBase>();

/**
 * This userobject collect values of a variable across an interface for each QP and compute a
 * scalar. The computed scalar value depends on the given parameter _interface_value_type\
 * _interface_value_type (see IntervafeValueTools).
 */
class InterfaceQpValueUserObjectBase : public InterfaceValueUserObject
{
public:
  InterfaceQpValueUserObjectBase(const InputParameters & parameters);
  virtual ~InterfaceQpValueUserObjectBase();

  virtual void initialize();
  virtual void execute() = 0;
  virtual void finalize() { return; };
  virtual void threadJoin(const UserObject & /*uo*/) { return; };

  virtual Real getQpValue(dof_id_type elem, unsigned int side, unsigned int qp) const;
  virtual Real getQpValueForLD(dof_id_type, unsigned int /*qp*/) const
  {
    mooseError("this function need to be overriden");
    return 0;
  };

protected:
  const bool _use_old_value;
  /// this map is used to store QP data.
  std::map<std::pair<dof_id_type, unsigned int>, std::vector<std::vector<Real>>> _map_values;
  // const VariableValue & _u;
  // const VariableValue & _u_neighbor;
  unsigned int _n_init = 0;
};
