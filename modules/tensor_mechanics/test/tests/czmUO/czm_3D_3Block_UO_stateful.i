[Mesh]
  file = coh3D_3Blocks.e
  parallel_type = REPLICATED
[]

[MeshModifiers]
  [./breakmesh]
    type = BreakMeshByBlock
    # split_interface = false
  [../]

  [./bottom_block_1]
    type = SideSetsAroundSubdomain
    depends_on = 'breakmesh'
    block = '1'
    new_boundary = 'bottom_1'
    normal = '0 0 -1'
  [../]
  [./top_block_2]
    type = SideSetsAroundSubdomain
    depends_on = 'breakmesh'
    block = '2'
    new_boundary = 'top_2'
    normal = '0 0 1'
  [../]
  [./top_block_3]
    type = SideSetsAroundSubdomain
    depends_on = 'breakmesh'
    block = '3'
    new_boundary = 'top_3'
    normal = '0 0 1'
  [../]

[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Modules/TensorMechanics/Master]
  [./all]
    strain = FINITE
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz stress_yz stress_xz stress_xy'
  [../]
[]


[BCs]
  [./bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = bottom_1
    value = 0.0
  [../]
  [./bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom_1
    value = 0.0
  [../]
  [./bottom_z]
    type = DirichletBC
    variable = disp_z
    boundary = bottom_1
    value = 0.0
  [../]
  [./top2_x]
    type = DirichletBC
    variable = disp_x
    boundary = top_2
    value = 0.0
  [../]
  [./top2_y]
    type = DirichletBC
    variable = disp_y
    boundary = top_2
    value = 0.0
  [../]
  [./top2_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = top_2
    function = 1*t
  [../]
  [./top3_x]
    type = DirichletBC
    variable = disp_x
    boundary = top_3
    value = 0.0
  [../]
  [./top3_y]
    type = DirichletBC
    variable = disp_y
    boundary = top_3
    value = 0.0
  [../]
  [./top3_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = top_3
    function = 1*t
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
    type = DispJumpAndNormalsUO_QP
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    boundary = 'interface'
    execute_on = 'initial NONLINEAR LINEAR timestep_end'
  [../]
  [./cohesive_law_exponential]
    type = CZMLaw3DC
    MaxAllowableTraction = '100 70'
    DeltaU0 = '1 0.7'
    displacement_jump_mp_name = 'displacement_jump_local'
    boundary = 'interface'
  [../]
[]

[Materials]
  [./Elasticity_tensor]
    type = ComputeElasticityTensor
    block = '1 2 3'
    fill_method = symmetric_isotropic
    C_ijkl = '0.3 0.5e8'
  [../]
  [./stress]
    type = ComputeFiniteStrainElasticStress
    block = '1 2 3'
  [../]
  [./gap]
    type = CZMUOBasedMaterial
    is_interface_material = true
    boundary = 'interface'
    displacement_jump_UO = 'displacement_jump_uo'
    traction_separation_UO = 'cohesive_law_exponential'
  [../]
  [./normal_MAT]
    type = CZMNormals
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
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  # petsc_options_value = 'hypre     boomerang'
  solve_type = NEWTON
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-6
  nl_max_its = 5
  l_tol = 1e-10
  l_max_its = 50
  start_time = 0.0
  dt = 0.2
  end_time = 5
  dtmin = 0.2
  line_search = none
  # num_steps = 1
[]
[Outputs]
  [./out]
    type = Exodus
  [../]
[]
[Postprocessors]
  [./sxx_3G]
    type = ElementAverageValue
    variable = stress_xx
    execute_on = 'initial timestep_end'
    block = 3
  [../]
  [./syy_3G]
    type = ElementAverageValue
    variable = stress_yy
    execute_on = 'initial timestep_end'
    block = 3
  [../]
  [./szz_3G]
    type = ElementAverageValue
    variable = stress_zz
    execute_on = 'initial timestep_end'
    block = 3
  [../]
  [./syz_3G]
    type = ElementAverageValue
    variable = stress_yz
    execute_on = 'initial timestep_end'
    block = 3
  [../]
  [./sxz_3G]
    type = ElementAverageValue
    variable = stress_xz
    execute_on = 'initial timestep_end'
    block = 3
  [../]
  [./sxy_3G]
    type = ElementAverageValue
    variable = stress_xy
    execute_on = 'initial timestep_end'
    block = 3
  [../]
  [./disp_3Z]
    type = ElementAverageValue
    variable = disp_z
    execute_on = 'initial timestep_end'
    block = 3
  [../]
[]
