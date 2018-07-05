/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef BREAKMESHBYBLOCKMANUAL_3BLOCKS_H
#define BREAKMESHBYBLOCKMANUAL_3BLOCKS_H

#include "BreakMeshByBlockManualBase.h"

// forward declaration
class BreakMeshByBlockManual_3Blocks;

template <>
InputParameters validParams<BreakMeshByBlockManual_3Blocks>();

class BreakMeshByBlockManual_3Blocks : public BreakMeshByBlockManualBase
{
public:
  BreakMeshByBlockManual_3Blocks(const InputParameters & parameters);

  virtual void modify() override;

private:
  void updateElements();
  void addInterfaceBoundary();
};

#endif // BREAKMESHBYBLOCKMANUAL_3BLOCKS_H
