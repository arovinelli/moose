[Mesh]
  [./msh]
    type = FileMeshGenerator
    file = coh3D_3Blocks.e
  []
  [./breakmesh]
    input = msh
    type = BreakMeshByBlockGenerator
  [../]
  [./bottom_block_1]
    input = breakmesh
    type = SideSetsAroundSubdomainGenerator
    block = '1'
    new_boundary = 'bottom_1'
    normal = '0 0 -1'
  [../]
  [./top_block_2]
    input = bottom_block_1
    type = SideSetsAroundSubdomainGenerator
    block = '2'
    new_boundary = 'top_2'
    normal = '0 0 1'
  [../]
  [./top_block_3]
    input = top_block_2
    type = SideSetsAroundSubdomainGenerator
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
    strain = SMALL
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
    type = FunctionDirichletBC
    variable = disp_x
    boundary = top_2
    function = 2*t
  [../]
  [./top2_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top_2
    function = 1*t
  [../]
  [./top2_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = top_2
    function = 3*t
  [../]
  [./top3_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = top_3
    function = 2*t
  [../]
  [./top3_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top_3
    function = 1*t
  [../]
  [./top3_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = top_3
    function = 3*t
  [../]
[]

[CohesiveZoneModel]
  boundary = 'interface'
  displacements = 'disp_x disp_y disp_z'
[]

[Materials]
  [./Elasticity_tensor]
    type = ComputeElasticityTensor
    block = '1 2 3'
    fill_method = symmetric_isotropic
    C_ijkl = '0.3 0.5e8'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = '1 2 3'
  [../]
  [./czm_3dc]
    type = CZM3DCLaw
    boundary = 'interface'
    normal_gap_at_maximum_normal_traction = 1
    tangential_gap_at_maximum_shear_traction = 0.5
    maximum_normal_traction = 100
    maximum_shear_traction = 70
    displacements = 'disp_x disp_y disp_z'
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
  dt = 0.1
  end_time = 1
  dtmin = 0.1
  line_search = none
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
  [./disp_top3_x]
    type = SideAverageValue
    variable = disp_x
    execute_on = 'initial timestep_end'
    boundary = 'top_3'
  [../]
  [./disp_top3_y]
    type = SideAverageValue
    variable = disp_y
    execute_on = 'initial timestep_end'
    boundary = 'top_3'
  [../]
  [./disp_top3_z]
    type = SideAverageValue
    variable = disp_z
    execute_on = 'initial timestep_end'
    boundary = 'top_3'
  [../]
[]
