[Mesh]
  file = testTri30.exo
  parallel_type = REPLICATED
[]

[MeshModifiers]
  [./breakmesh]
    type = BreakMeshByBlock
    # split_interface = false
  [../]
  [./add_side_sets]
     type = SideSetsFromNormals
     normals = '0  -1  0
                0  1  0
                -1 0  0
                1  0  0
                0 0 -1
                0 0  1'
     fixed_normal = true
     new_boundary = 'bottom top left right rear front'
     depends_on = breakmesh
   [../]
   # [./left_n1]
   #   type = AddExtraNodeset
   #   new_boundary = 'left_n1'
   #   nodes = '2 69 114 191 246'
   #   depends_on = breakmesh
   # [../]
   # [./left_n2]
   #   type = AddExtraNodeset
   #   new_boundary = 'left_n2'
   #   nodes = '7 37 25 2 45'
   #   depends_on = breakmesh
   # [../]
   # [./right_n1]
   #   type = AddExtraNodeset
   #   new_boundary = 'right_n1'
   #   nodes = '225 288 310 383'
   #   depends_on = breakmesh
   # [../]
   # [./right_n2]
   #   type = AddExtraNodeset
   #   new_boundary = 'right_n2'
   #   nodes = '228 276 296 533'
   #   depends_on = breakmesh
   # [../]
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Modules/TensorMechanics/Master]
  [./all]
    strain = SMALL
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz stress_yz stress_xz stress_xy'
  [../]
[]


[BCs]
  # [./left_x]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = left
  #   value = 0.0
  # [../]
  # [./left_y]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = left
  #   value = 0.0
  # [../]
  #
  # [./right_x]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = right
  #   value = 0.0
  # [../]
  # [./right_y]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = right
  #   value = 0.0
  # [../]
  #
  # [./top_x]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = top
  #   value = 0.0
  # [../]
  # [./top_y]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = top
  #   value = 0.0
  # [../]
  #
  # [./bottom_x]
  #   type = DirichletBC
  #   variable = disp_x
  #   boundary = bottom
  #   value = 0.0
  # [../]
  # [./bottom_y]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = bottom
  #   value = 0.0
  # [../]

  [./front_x]
    type = DirichletBC
    variable = disp_x
    boundary = front
    value = 0.0
  [../]
  [./front_y]
    type = DirichletBC
    variable = disp_y
    boundary = front
    value = 0.0
  [../]
  [./front_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = front
    function = loadUnloadFunction
  [../]

  [./rear_x]
    type = DirichletBC
    variable = disp_x
    boundary = front
    value = 0.0
  [../]
  [./rear_y]
    type = DirichletBC
    variable = disp_y
    boundary = front
    value = 0.0
  [../]
  [./rear_z]
    type = DirichletBC
    variable = disp_z
    boundary = rear
    value = 0.0
  [../]
[]
# [ICs]
#   [./ic_x]
#     type = RandomIC
#     variable = disp_x
#     boundary = 'interface'
#     seed = 0
#     min = -0.1
#     max =  0.00
#     legacy_generator = false
#   [../]
#   [./ic_y]
#     type = RandomIC
#     variable = disp_y
#     boundary = 'interface'
#     min = -0.1
#     max =  0.00
#     seed = 0
#     legacy_generator = false
#   [../]
#   [./ic_z]
#     type = RandomIC
#     variable = disp_z
#     boundary = 'interface'
#     min = -0.1
#     max =  0.00
#     seed = 0
#     legacy_generator = false
#   [../]
# []
[Functions]
  [./loadUnloadFunction]
    type = PiecewiseLinear
    x = '0 4      8  14     20  21     32    42   67   92 142'
    y = '0 0.00002  0  0.00012  0   -0.0002   0.002   0    0.005   0   0.1'
    # x = '0 4    '
    # y = '0 -0.002'
  [../]
[]
[InterfaceKernels]
  [./interface_x]
    type = CZMInterfaceKernel
    variable = disp_x
    neighbor_var = disp_x
    disp_1 = disp_y
    disp_1_neighbor = disp_y
    disp_2 = disp_z
    disp_2_neighbor = disp_z
    disp_index = 0
    boundary = 'interface'
  [../]
  [./interface_y]
    type = CZMInterfaceKernel
    variable = disp_y
    neighbor_var = disp_y
    disp_1 = disp_x
    disp_1_neighbor = disp_x
    disp_2 = disp_z
    disp_2_neighbor = disp_z
    disp_index = 1
    boundary = 'interface'
  [../]
  [./interface_z]
    type = CZMInterfaceKernel
    variable = disp_z
    neighbor_var = disp_z
    disp_1 = disp_x
    disp_1_neighbor = disp_x
    disp_2 = disp_y
    disp_2_neighbor = disp_y
    disp_index = 2
    boundary = 'interface'
  [../]
[]
[UserObjects]
  [./displacement_jump_uo]
    type = DispJumpUO_QP
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    boundary = 'interface'
    execute_on = 'initial LINEAR NONLINEAR timestep_end'
  [../]

  [./cohesive_law_exponential]
    type = CZMLawExponential
    displacement_jump_peak = 0.001
    traction_peak = 150
    displacement_jump_mp_name = 'displacement_jump_local'
    boundary = 'interface'
    compression_multiplier = 1e3
  [../]
  # [./cohesive_law_unload_linear]
  #   type = CZMUnloadLinear
  #   displacement_jump_mp_name = 'displacement_jump_local'
  #   boundary = 'interface'
  # [../]
  # [./cohesive_law_copenetration]
  #   type = CZMCopenetrationPenalty
  #   displacement_jump_mp_name = 'displacement_jump_local'
  #   boundary = 'interface'
  #   copenetration_penalty_stiffness = 2e7
  # [../]
[]



[Materials]
  [./Elasticity_tensor]
    type = ComputeElasticityTensor
    fill_method = symmetric_isotropic
    C_ijkl = '0.3 2e5'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
  [./gap]
    type = CZMUOBasedMaterial
    is_interface_material = true
    boundary = 'interface'
    displacement_jump_UO = 'displacement_jump_uo'
    traction_separation_UO = 'cohesive_law_exponential'
    # unload_traction_separation_UO = 'cohesive_law_unload_linear'
    # coopenetration_penalty_UO = 'cohesive_law_copenetration'
    # coopenetration_penalty = 100
  [../]
[]
 [Preconditioning]
   [./SMP]
     type = SMP
     full = true
   [../]
 []
[Executioner]
  # Preconditisoned JFNK (default)
  type = Transient
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu        superlu_dist'
  # petsc_options_value = 'hypre     boomerang'
  solve_type = NEWTON
  nl_abs_tol = 1e-5
  nl_rel_tol = 1e-6
  nl_max_its = 20
  l_tol = 1e-15
  # l_max_its = 10
  start_time = 0.0
  dt = 1
  # dtmin = 0.01
  end_time = 142
  line_search = none
[]
[Outputs]
  [./out]
    type = Exodus
  [../]
[]
[Postprocessors]
  [./sxx]
    type = ElementAverageValue
    variable = stress_xx
    execute_on = 'initial timestep_end'
  [../]
  [./syy]
    type = ElementAverageValue
    variable = stress_yy
    execute_on = 'initial timestep_end'
  [../]
  [./szz]
    type = ElementAverageValue
    variable = stress_zz
    execute_on = 'initial timestep_end'
  [../]
  [./syz]
    type = ElementAverageValue
    variable = stress_yz
    execute_on = 'initial timestep_end'
  [../]
  [./sxz]
    type = ElementAverageValue
    variable = stress_xz
    execute_on = 'initial timestep_end'
  [../]
  [./sxy]
    type = ElementAverageValue
    variable = stress_xy
    execute_on = 'initial timestep_end'
  [../]
  [./disp_rightZ]
    type = SideAverageValue
    variable = disp_z
    execute_on = 'initial timestep_end'
    boundary = 'front'
  [../]
[]
