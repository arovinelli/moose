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

#ifndef BREAKMESHBYBLOCKMANUAL_2DJUNCTION_H
#define BREAKMESHBYBLOCKMANUAL_2DJUNCTION_H

#include "BreakMeshByBlockManualBase.h"

// forward declaration
class BreakMeshByBlockManual_2DJunction;

template <>
InputParameters validParams<BreakMeshByBlockManual_2DJunction>();

class BreakMeshByBlockManual_2DJunction : public BreakMeshByBlockManualBase
{
public:
  BreakMeshByBlockManual_2DJunction(const InputParameters & parameters);

  virtual void modify() override;

private:
  void updateElements();
  void addInterfaceBoundary();
};

#endif // BREAKMESHBYBLOCKMANUAL_2DJUNCTION_H
