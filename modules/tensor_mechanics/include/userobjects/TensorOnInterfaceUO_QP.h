//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef TENSORONINTERFACEUO_QP
#define TENSORONINTERFACEUO_QP

#include "InterfaceUserObject.h"

class TensorOnInterfaceUO_QP;

template <>
InputParameters validParams<TensorOnInterfaceUO_QP>();

/**
 *
 */
class TensorOnInterfaceUO_QP : public InterfaceUserObject
{
public:
  TensorOnInterfaceUO_QP(const InputParameters & parameters);
  virtual ~TensorOnInterfaceUO_QP();

  virtual void initialize();
  virtual void execute();
  virtual void finalize() { return; };
  virtual void threadJoin(const UserObject & /*uo*/) { return; };

  RankTwoTensor getTensorMaster(dof_id_type elem, unsigned int side, unsigned int qp) const;
  RankTwoTensor getTensorSlave(dof_id_type elem, unsigned int side, unsigned int qp) const;

protected:
  /// this map is used for storing data at QPs.
  /// keys<element_id, side_id>, values<vector (1 elem_per_QP)<vector<Real (_mean_mat_prop, _var_jump)>>>
  std::map<std::pair<dof_id_type, unsigned int>, std::vector<std::vector<RankTwoTensor>>>
      _map_values;
  const MaterialProperty<RankTwoTensor> & _tensor_M;
  const MaterialProperty<RankTwoTensor> & _tensor_S;
};

#endif // TENSORONINTERFACEUO_QP
