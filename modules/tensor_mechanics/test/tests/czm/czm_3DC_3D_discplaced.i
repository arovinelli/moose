[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [./msh]
    type = GeneratedMeshGenerator
    nx = 1
    ny = 1
    nz = 4
    dim = 3
  []
  [./subdomain_1]
    type = SubdomainBoundingBoxGenerator
    input = msh
    bottom_left = '0 0 0'
    block_id = 1
    top_right = '1 1 0.5'
  []
  [./subdomain_2]
    type = SubdomainBoundingBoxGenerator
    input = subdomain_1
    bottom_left = '0 0 0.5'
    block_id = 2
    top_right = '1 1 1'
  []
  [./breakmesh]
    input = subdomain_2
    type = BreakMeshByBlockGenerator
  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    strain = SMALL
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz stress_yz stress_xz stress_xy'
  [../]
[]

[InterfaceKernels]

  [./czmz]
    type = CZMInterfaceKernel
    boundary = 'interface'
    displacements = 'disp_x disp_y disp_z'
    use_displaced_mesh = true
    component = 2
    save_in = 'R_z R_neigh_z'
    save_in_var_side = 'm s'
    variable = disp_z
    neighbor_var = disp_z
  [../]
[]

[BCs]
  [./left_x]
    type = DirichletBC
    variable = disp_x
    preset = true
    boundary = left
    value = 0.0
  [../]
  [./left_y]
    type = DirichletBC
    variable = disp_y
    preset = true
    boundary = left
    value = 0.0
  [../]
  [./left_z]
    type = DirichletBC
    variable = disp_z
    preset = true
    boundary = back
    value = 0.0
  [../]
  [./right_z]
    type = FunctionDirichletBC
    variable = disp_z
    preset = true
    boundary = front
    function = 1*t
  [../]
[]
[AuxVariables]
[./R_x]
  family = LAGRANGE
  order = FIRST
[]
[./R_y]
  family = LAGRANGE
  order = FIRST
[]
[./R_z]
  family = LAGRANGE
  order = FIRST
[]
[./R_neigh_x]
  family = LAGRANGE
  order = FIRST
[]
[./R_neigh_y]
  family = LAGRANGE
  order = FIRST
[]
[./R_neigh_z]
  family = LAGRANGE
  order = FIRST
[]
[]

[Materials]
  [./Elasticity_tensor]
    type = ComputeElasticityTensor
    block = '1 2'
    fill_method = symmetric_isotropic
    C_ijkl = '0.3 0.5e8'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = '1 2'
  [../]
  [./czm_3dc]
    type = SalehaniIrani3DCTraction
    boundary = 'interface'
    normal_gap_at_maximum_normal_traction = 1
    tangential_gap_at_maximum_shear_traction = 0.5
    maximum_normal_traction = 100
    maximum_shear_traction = 70
    displacements = 'disp_x disp_y disp_z'
    use_displaced_mesh = true
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
  nl_rel_tol = 1e-6
  nl_max_its = 5
  l_tol = 1e-10
  l_max_its = 50
  start_time = 0.0
  dt = 0.2
  end_time = 1
  dtmin = 0.2
  line_search = none
[]

[Outputs]
  [./out]
    type = Exodus
  [../]
[]

[Postprocessors]
  [./sxx]
    type = SideAverageValue
    variable = stress_xx
    execute_on = 'INITIAL TIMESTEP_END'
    boundary = 'interface'
  [../]
  [./syy]
    type = SideAverageValue
    variable = stress_yy
    execute_on = 'INITIAL TIMESTEP_END'
    boundary = 'interface'
  [../]
  [./szz]
    type = SideAverageValue
    variable = stress_zz
    execute_on = 'INITIAL TIMESTEP_END'
    boundary = 'interface'
  [../]
  [./syz]
    type = SideAverageValue
    variable = stress_yz
    execute_on = 'INITIAL TIMESTEP_END'
    boundary = 'interface'
  [../]
  [./sxz]
    type = SideAverageValue
    variable = stress_xz
    execute_on = 'INITIAL TIMESTEP_END'
    boundary = 'interface'
  [../]
  [./sxy]
    type = SideAverageValue
    variable = stress_xy
    execute_on = 'INITIAL TIMESTEP_END'
    boundary = 'interface'
  [../]
  [./disp_x]
    type = SideAverageValue
    variable = disp_x
    execute_on = 'INITIAL TIMESTEP_END'
    boundary = 'right'
  [../]
  [./disp_y]
    type = SideAverageValue
    variable = disp_y
    execute_on = 'INITIAL TIMESTEP_END'
    boundary = 'right'
  [../]
  [./disp_z]
    type = SideAverageValue
    variable = disp_z
    execute_on = 'INITIAL TIMESTEP_END'
    boundary = 'right'
  [../]
[]
