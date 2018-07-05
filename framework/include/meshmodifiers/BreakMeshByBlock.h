//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef BREAKMESHBYBLOCK_H
#define BREAKMESHBYBLOCK_H

#include "BreakMeshByBlockBase.h"

class BreakMeshByBlock;

template <>
InputParameters validParams<BreakMeshByBlock>();

class BreakMeshByBlock : public BreakMeshByBlockBase
{
public:
  BreakMeshByBlock(const InputParameters & parameters);

  virtual void modify() override;

private:
  /// loop over the node and find their support
  void buildNodeSupport();

  /// find nodes duplicity (i.e. how many times a node need to be duplicated)
  void buildInterfacialNodes();

  /// duplicate nodes with duplicity greater than 1
  void duplicateNodes();

  /// assign duplicated nodes to the approriate element
  void tearElements();

  /// generate the new boundary interface
  void addInterfaceBoundary();

  std::map<dof_id_type, std::vector<dof_id_type>> _node_support;
  std::set<std::pair<subdomain_id_type, subdomain_id_type>> _boundary_pairs;
  std::map<std::pair<subdomain_id_type, subdomain_id_type>,
           std::set<std::pair<dof_id_type, unsigned int>>>
      _boundary_sides;
  std::map<dof_id_type, unsigned int> _node_duplicity;
  std::map<dof_id_type, std::vector<subdomain_id_type>> _duplicated_node_materials;
  std::map<dof_id_type, std::vector<dof_id_type>> _duplicated_node;
};

#endif /* BREAKMESHBYBLOCK_H */
