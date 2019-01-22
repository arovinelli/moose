//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SCALARBULKMPACROSSINTERFACE_QP
#define SCALARBULKMPACROSSINTERFACE_QP

#include "InterfaceUserObject.h"

class ScalarBulkMPAcrossInterface_QP;

template <>
InputParameters validParams<ScalarBulkMPAcrossInterface_QP>();

/**
 *
 */
class ScalarBulkMPAcrossInterface_QP : public InterfaceUserObject
{
public:
  ScalarBulkMPAcrossInterface_QP(const InputParameters & parameters);
  virtual ~ScalarBulkMPAcrossInterface_QP();

  virtual void initialize();
  virtual void execute();
  virtual void finalize() { return; };
  virtual void threadJoin(const UserObject & /*uo*/) { return; };

  Real getValueAverage(dof_id_type elem, unsigned int side, unsigned int qp) const;
  Real getValueMaster(dof_id_type elem, unsigned int side, unsigned int qp) const;
  Real getValueSlave(dof_id_type elem, unsigned int side, unsigned int qp) const;
  unsigned int getScalarVarID() const { return _scalar_var_id; };

protected:
  /// this map is used for storing data at QPs.
  /// keys<element_id, side_id>, values<vector (1 elem_per_QP)<vector<Real
  /// (_mean_mat_prop, _var_jump)>>>
  std::map<std::pair<dof_id_type, unsigned int>, std::vector<std::vector<Real>>> _map_values;
  const VariableValue & _scalar_M;
  const VariableValue & _scalar_S;
  const unsigned int _scalar_var_id;
};

#endif // SCALARBULKMPACROSSINTERFACE_QP
