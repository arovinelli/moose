[Mesh]
  [./msh]
    type = GeneratedMeshGenerator
    nx =2
    dim = 1
  []
  [./subdomain_1]
    type = SubdomainBoundingBoxGenerator
    input = msh
    bottom_left = '0 0 0'
    block_id = 1
    top_right = '0.5 0 0'
  []
  [./subdomain_2]
    type = SubdomainBoundingBoxGenerator
    input = subdomain_1
    bottom_left = '0.5 0 0'
    block_id = 2
    top_right = '1 0 0'
  []
  [./breakmesh]
    input = subdomain_2
    type = BreakMeshByBlockGenerator
  [../]
[]
[GlobalParams]
  displacements = 'disp_x'
[]

[Modules/TensorMechanics/Master]
  [./all]
    strain = SMALL
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz stress_yz stress_xz stress_xy'
  [../]
[]

[BCs]
  [./left_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./right_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = 1*t
  [../]
[]

[InterfaceKernels]
  [./interface_x]
    type = CZMInterfaceKernel
    variable = disp_x
    neighbor_var = disp_x
    displacements = 'disp_x'
    component = 0
    boundary = 'interface'
  [../]
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
    type = CZM3DCLaw
    boundary = 'interface'
    normal_gap_at_maximum_normal_traction = 1
    tangential_gap_at_maximum_shear_traction = 0.5
    maximum_normal_traction = 100
    maximum_shear_traction = 70
    displacements = 'disp_x'
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
  end_time = 5
  dtmin = 0.2
  line_search = none
[]

[Outputs]
  [./out]
    type = Exodus
  [../]
[]

[Postprocessors]
  [./sxx_2G]
    type = ElementAverageValue
    variable = stress_xx
    execute_on = 'initial timestep_end'
    block = 2
  [../]
  [./syy_2G]
    type = ElementAverageValue
    variable = stress_yy
    execute_on = 'initial timestep_end'
    block = 2
  [../]
  [./szz_2G]
    type = ElementAverageValue
    variable = stress_zz
    execute_on = 'initial timestep_end'
    block = 2
  [../]
  [./syz_2G]
    type = ElementAverageValue
    variable = stress_yz
    execute_on = 'initial timestep_end'
    block = 2
  [../]
  [./sxz_2G]
    type = ElementAverageValue
    variable = stress_xz
    execute_on = 'initial timestep_end'
    block = 2
  [../]
  [./sxy_2G]
    type = ElementAverageValue
    variable = stress_xy
    execute_on = 'initial timestep_end'
    block = 2
  [../]
  [./disp_right_x]
    type = SideAverageValue
    variable = disp_x
    execute_on = 'initial timestep_end'
    boundary = 'right'
  [../]
[]