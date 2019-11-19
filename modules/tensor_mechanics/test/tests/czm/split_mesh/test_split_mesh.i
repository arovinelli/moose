[Mesh]
  [./msh]
    type = FileMeshGenerator
    file = testTri10.e
  []

  [./breakmesh]
    input = msh
    type = BreakMeshByBlockGenerator
  [../]

  [./add_side_sets]
    input = breakmesh
    type = SideSetsFromNormalsGenerator
    normals = '0 -1  0
               0  1  0
               -1 0  0
               1  0  0
               0  0 -1
               0  0  1'
    fixed_normal = true
    new_boundary = 'y0 y1 x0 x1 z0 z1'
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

[Modules/TensorMechanics/CohesiveZoneMaster]
  [./czm1]
    boundary = 'interface'
    displacements = 'disp_x disp_y disp_z'
  [../]
[]

[Functions]
  [./applied_load_x]
    type = PiecewiseLinear
    x = '0 0.1 1e6'
    y = '0 0 0'
  [../]
  [./applied_load_y]
    type = PiecewiseLinear
    x = '0 0.1 1e6'
    y = '0 0 0'
  [../]
  [./applied_load_z]
    type = PiecewiseLinear
    x = '0 0.1 1e6'
    y = '0 180 180'
  [../]
[]

[BCs]
  [./x0]
    type = DirichletBC
    variable = disp_x
    boundary = x0
    value = 0.0
  [../]
  [./y0]
    type = DirichletBC
    variable = disp_y
    boundary = y0
    value = 0.0
  [../]
  [./z0]
    type = DirichletBC
    variable = disp_z
    boundary = z0
    value = 0.0
  [../]
  [./x1]
    type = FunctionNeumannBC
    boundary = x1
    function = applied_load_x
    variable = disp_x
  [../]
  [./y1]
    type = FunctionNeumannBC
    boundary = y1
    function = applied_load_y
    variable = disp_y
  [../]
  [./z1]
    type = FunctionNeumannBC
    boundary = z1
    function = applied_load_z
    variable = disp_z
  [../]
[]

# Constraint System
[Constraints]
  [./x1]
    type = EqualValueBoundaryConstraint
    variable = disp_x
    slave =  x1    # boundary
    penalty = 1e6
  [../]
  [./y1]
    type = EqualValueBoundaryConstraint
    variable = disp_y
    slave = y1    # boundary
    penalty = 1e6
  [../]
  [./z1]
    type = EqualValueBoundaryConstraint
    variable = disp_z
    slave = z1    # boundary
    penalty = 1e6
  [../]
[]


[Materials]
  [./Elasticity_tensor]
    type = ComputeElasticityTensor
    fill_method = symmetric_isotropic
    C_ijkl = '0.3 0.5e8'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
  [./czm_3dc]
    type = SalehaniIrani3DCTraction
    boundary = 'interface'
    normal_gap_at_maximum_normal_traction = 1
    tangential_gap_at_maximum_shear_traction = 0.5
    maximum_normal_traction = 1000
    maximum_shear_traction = 700
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
  dt = 0.2
  end_time = 0.2
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
    execute_on = 'initial timestep_end'
    boundary = 'interface'
  [../]
  [./syy]
    type = SideAverageValue
    variable = stress_yy
    execute_on = 'initial timestep_end'
    boundary = 'interface'
  [../]
  [./szz]
    type = SideAverageValue
    variable = stress_zz
    execute_on = 'initial timestep_end'
    boundary = 'interface'
  [../]
  [./syz]
    type = SideAverageValue
    variable = stress_yz
    execute_on = 'initial timestep_end'
    boundary = 'interface'
  [../]
  [./sxz]
    type = SideAverageValue
    variable = stress_xz
    execute_on = 'initial timestep_end'
    boundary = 'interface'
  [../]
  [./sxy]
    type = SideAverageValue
    variable = stress_xy
    execute_on = 'initial timestep_end'
    boundary = 'interface'
  [../]
[]
