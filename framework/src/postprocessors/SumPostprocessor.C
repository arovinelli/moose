//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SumPostprocessor.h"

registerMooseObject("MooseApp", SumPostprocessor);

template <>
InputParameters
validParams<SumPostprocessor>()
{
  InputParameters params = validParams<GeneralPostprocessor>();
  params.addRequiredParam<PostprocessorName>("value1", "First value");
  params.addRequiredParam<PostprocessorName>("value2", "Second value");

  return params;
}

SumPostprocessor::SumPostprocessor(const InputParameters & parameters)
  : GeneralPostprocessor(parameters),
    _value1(getPostprocessorValue("value1")),
    _value2(getPostprocessorValue("value2"))
{
}

void
SumPostprocessor::initialize()
{
}

void
SumPostprocessor::execute()
{
}

PostprocessorValue
SumPostprocessor::getValue()
{
  return _value1 + _value2;
}
