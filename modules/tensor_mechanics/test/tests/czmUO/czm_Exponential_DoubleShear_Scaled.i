[Mesh]
  file = 3D_3Block_3x3.e
  parallel_type = REPLICATED
[]

[MeshModifiers]
  [./scale]
    type = Transform
    transform = SCALE
    vector_value = '2 2 1'
  []
  [./breakmesh]
    type = BreakMeshByBlock
    depends_on = scale
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

[Functions]
  [./loadUnloadFunction]
    type = PiecewiseLinear
    x = '0 10 20 220'
    y = '0 1.41421356237  0  -28.2842712475'
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
    type = FunctionDirichletBC
    variable = disp_x
    boundary = top_2
    function = loadUnloadFunction
  [../]
  [./top2_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top_2
    function = loadUnloadFunction
  [../]
  [./top2_z]
    type = DirichletBC
    variable = disp_z
    boundary = top_2
    value = 0.0
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

[AuxVariables]
  [./Tn_interface]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Tt_interface]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]
[AuxKernels]
  [./Tn]
    type = MaterialRealVectorValueAux
    component = 0
    property = traction_local
    variable = Tn_interface
    boundary = interface
  []
  [./Tt]
    type = MaterialRealAux
    property = shear_traction
    variable = Tt_interface
    boundary = interface
  []
[]


[UserObjects]
  [./displacement_jump_uo]
    type = DispJumpUO_QP
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    boundary = 'interface'
    execute_on = 'initial LINEAR timestep_end'
  [../]
  [./cohesive_law_exponential]
    type = CZMLawExponential
    displacement_jump_peak = 1
    traction_peak = 100
    displacement_jump_mp_name = 'displacement_jump_local'
    boundary = 'interface'
    compression_multiplier = 1e3
    beta = 0.5
  [../]
[]

[Materials]
  [./Elasticity_tensor]
    type = ComputeElasticityTensor
    block = '1 2'
    fill_method = symmetric_isotropic
    C_ijkl = '0.3 200e3'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = '1 2'
  [../]
  [./gap]
    type = CZMUOBasedMaterial
    is_interface_material = true
    boundary = 'interface'
    displacement_jump_UO = 'displacement_jump_uo'
    traction_separation_UO = 'cohesive_law_exponential'
    compute_shear_traction = true
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
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  solve_type = NEWTON
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-8
  nl_max_its = 5

  l_max_its = 50
  start_time = 0.0
  dt = 10
  end_time = 220

  line_search = none
[]

[Outputs]
  [./out]
    type = Exodus
  [../]
[]

[Postprocessors]
  [./sxx_2G]
    type = SideAverageValue
    variable = stress_xx
    execute_on = 'initial timestep_end'
    boundary = 'top_2'
  [../]
  [./syy_2G]
    type = SideAverageValue
    variable = stress_yy
    execute_on = 'initial timestep_end'
    boundary = 'top_2'
  [../]
  [./szz_2G]
    type = SideAverageValue
    variable = stress_zz
    execute_on = 'initial timestep_end'
    boundary = 'top_2'
  [../]
  [./syz_2G]
    type = SideAverageValue
    variable = stress_yz
    execute_on = 'initial timestep_end'
    boundary = 'top_2'
  [../]
  [./sxz_2G]
    type = SideAverageValue
    variable = stress_xz
    execute_on = 'initial timestep_end'
    boundary = 'top_2'
  [../]
  [./sxy_2G]
    type = SideAverageValue
    variable = stress_xy
    execute_on = 'initial timestep_end'
    boundary = 'top_2'
  [../]


  [./disp_top2_z]
    type = SideAverageValue
    variable = disp_z
    execute_on = 'initial timestep_end'
    boundary = 'top_2'
  [../]
  [./disp_top2_y]
    type = SideAverageValue
    variable = disp_y
    execute_on = 'initial timestep_end'
    boundary = 'top_2'
  [../]
  [./disp_top2_x]
    type = SideAverageValue
    variable = disp_x
    execute_on = 'initial timestep_end'
    boundary = 'top_2'
  [../]
[]
