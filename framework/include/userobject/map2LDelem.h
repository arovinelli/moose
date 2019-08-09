//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef map2LDelem_h
#define map2LDelem_h

#include "InterfaceUserObject.h"

class map2LDelem;

template <>
InputParameters validParams<map2LDelem>();

/**
 *
 */
class map2LDelem : public InterfaceUserObject
{
public:
  map2LDelem(const InputParameters & parameters);
  virtual ~map2LDelem();

  virtual void initialize() override;
  virtual void execute() override{};
  virtual void finalize() override { return; };
  virtual void threadJoin(const UserObject & /*uo*/) override { return; };
  // Real getQpValueForLD(dof_id_type ld_elem, unsigned int qp) const;

  std::pair<dof_id_type, unsigned int> getLDNeighbor(dof_id_type ld_elem) const;

protected:
  const std::vector<SubdomainID> ld_block_ids;
  // mapping element side to lower element id
  std::map<dof_id_type, std::pair<dof_id_type, unsigned int>> _map_LD_with_elem_side;
};

#endif // map2LDelem_h
