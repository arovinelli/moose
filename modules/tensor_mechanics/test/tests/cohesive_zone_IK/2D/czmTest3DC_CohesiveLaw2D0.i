[Mesh]
  type = CohesiveZoneMeshSplit
  file = czm2Blocks.e
  displacements = 'disp_x disp_y '
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]




[AuxVariables]
  [./sxx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./syy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./sxy]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]

  [./sxx]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    variable = sxx
    block = '1 2'
  []
  [./syy]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    variable = syy
    block = '1 2'
  []
  [./sxy]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    variable = sxy
    block = '1 2'
  []

[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y '
  [../]
[]

[BCs]
  [./bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = 1
    value = 0.0
  [../]
  [./bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = 1
    value = 0.0
  [../]


  [./top_x]
    type = DirichletBC
    variable = disp_x
    boundary = 2
    value = 0.0
  [../]
  [./top_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 2
    function = 0.01*t
  [../]




[]

[InterfaceKernels]
  [./interface_x]
    type = DisplacementJumpBasedCohesiveInterfaceKernel
    variable = disp_x
    neighbor_var = disp_x
    disp_1 = disp_y
    disp_1_neighbor = disp_y
    disp_index = 0
    boundary = 100

  [../]
  [./interface_y]
    type = DisplacementJumpBasedCohesiveInterfaceKernel
    variable = disp_y
    neighbor_var = disp_y
    disp_1 = disp_x
    disp_1_neighbor = disp_x
    disp_index = 1
    boundary = 100
  [../]

[]

[UserObjects]
  # [./TractionSeparationTest]
  #   type = CohesiveLaw_3DC
  #   boundary = 100
  #   DeltaU0 = '1 0.5'
  #   MaxAllowableTraction = '1e2 5e1'
  # [../]
  [./DisplacementJump]
    type = DisplacementJumpCohesiveInterface
    # boundary = 100

    disp_x = disp_x
    disp_x_neighbor = disp_x
    disp_y = disp_y
    disp_y_neighbor = disp_y
    execute_on = 'LINEAR'
  [../]
[]

[Materials]
  [./Elasticity_tensor]
    type = ComputeElasticityTensor
    block = '1 2'
    fill_method = symmetric_isotropic
    C_ijkl = '0.3 0.5e8'
  [../]
  [./strain]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y '
    block = '1 2'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = '1 2'
  [../]
  [./gap]
    type = DisplacementJumpBasedCohesiveInterfaceMaterial
    boundary = 100
    # uo_TractionSeparationLaw = 'TractionSeparationTest'
    uo_CohesiveInterface = 'DisplacementJump'

  [../]
[]

 [Preconditioning]
   [./SMP]
     type = SMP
     full = true
   [../]
 []

[Executioner]
  type = Transient

  # petsc_options = '-snes_linesearch_monitor -snes_monitor -ksp_monitor_true_residual -snes_converged_reason -ksp_converged_reason nes_grid_sequence'
  petsc_options_iname = '-pc_type '
  petsc_options_value = 'lu '

  solve_type = NEWTON
  nl_abs_tol = 1e-6
  nl_rel_tol = 1e-8
  nl_max_its = 50
  l_tol = 1e-10
  # line_search = none
  l_max_its = 50
  start_time = 0.0
  dt = 100
  num_steps = 1
  # end_time = 1000

[]

[Outputs]
  [./out]
    type = Exodus
  [../]
[]

[Postprocessors]
  [./sxx_G]
    type = ElementAverageValue
    variable = sxx
    execute_on = 'initial timestep_end'
    block = 1
  [../]
  [./syy_G]
    type = ElementAverageValue
    variable = syy
    execute_on = 'initial timestep_end'
    block = 1
  [../]

  [./sxy_G]
    type = ElementAverageValue
    variable = sxy
    execute_on = 'initial timestep_end'
    block = 1
  [../]


[]
