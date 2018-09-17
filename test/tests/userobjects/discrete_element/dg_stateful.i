[Mesh]
  dim = 3
  file = cube.e
[]

[Variables]
  [./u]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./prop1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./A_0]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./A_1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./A_0_old]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./A_1_old]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./B_0]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./B_1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./B_2]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./heat]
    type = MatDiffusionTest
    variable = u
    prop_name = thermal_conductivity
    prop_state = 'old'                  # Use the "Old" value to compute conductivity
  [../]
  [./ie]
    type = TimeDerivative
    variable = u
  [../]
[]

[AuxKernels]
  [./prop1_output]
    type = MaterialRealAux
    variable = prop1
    property = thermal_conductivity
  [../]

  [./prop1_output_init]
    type = MaterialRealAux
    variable = prop1
    property = thermal_conductivity
    execute_on = initial
  [../]

  [./A_0]
    type = MaterialStdVectorAux
    variable = A_0
    property = A
    index = 0
    execute_on = timestep_end
  [../]
  [./A_1]
    type = MaterialStdVectorAux
    variable = A_1
    property = A
    index = 1
    execute_on = timestep_end
  [../]
  [./A_0_old]
    type = MaterialStdVectorAux
    variable = A_0_old
    property = A
    index = 0
    execute_on = timestep_begin

  [../]
  [./A_1_old]
    type = MaterialStdVectorAux
    variable = A_1_old
    property = A
    index = 1
    execute_on = timestep_begin
  [../]


  [./B_0]
    type = MaterialStdVectorAux
    variable = B_0
    property = B
    index = 0
    execute_on = 'timestep_end'
  [../]
  [./B_1]
    type = MaterialStdVectorAux
    variable = B_1
    property = B
    index = 1
    execute_on = 'timestep_end'
  [../]
  [./B_2]
    type = MaterialStdVectorAux
    variable = B_2
    property = B
    index = 2
    execute_on = 'timestep_end'
  [../]
[]

[BCs]
  [./bottom]
    type = DirichletBC
    variable = u
    boundary = 1
    value = 0.0
  [../]
  [./top]
    type = DirichletBC
    variable = u
    boundary = 2
    value = 1.0
  [../]
[]

[UserObjects]
  [test_discrete_UO]
    type = DiscreteElementWritingMP
    MP_names = 'A B'
    MP_sizes = '2 3'
    MP_stateful = '1 0'
    MP_initial_values = '1. 2.
                         3. 4. 5.'
  [../]
[../]

[Materials]
  [./stateful]
    type = StatefulTest
    prop_names = thermal_conductivity
    prop_values = 1.0
  [../]
  [./DEUOmaterial]
    type = DiscreteElementUOMaterial
    discrete_elem_uo = 'test_discrete_UO'
  [../]
[]

[Postprocessors]
  [./integral]
    type = ElementAverageValue
    variable = prop1
    execute_on = 'initial timestep_end'
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'
  l_max_its = 10
  start_time = 0.0
  num_steps = 5
  dt = .1
[]

[Outputs]
  file_base = out
  exodus = true
[]
