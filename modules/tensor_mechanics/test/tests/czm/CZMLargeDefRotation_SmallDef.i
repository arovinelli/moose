# Simple 3D test

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
      [./disp_x]
      [../]
      [./disp_y]
      [../]
      [./disp_z]
      [../]
[]

[Mesh]
  [./msh]
  type = FileMeshGenerator
  file = patch_mesh.e
  []
  [./transform]
  type = TransformGenerator
  input = msh
  transform = TRANSLATE
  vector_value = '-0.5 -0.5 -0.5'
  []
  [./split]
    type = BreakMeshByBlockGenerator
    input = transform
  []
  [./add_surfaces]
    type = SideSetsFromNormalsGenerator
    input = split
    normals = '0  0  1
               0  1  0
               1  0  0
               0  0 -1
               0 -1  0
              -1  0  0'
    fixed_normal = true
    new_boundary = 'front top right back bottom left'
  []
[]

[Modules]
  [./TensorMechanics]
    [./Master]
      [./all]
        strain = FINITE
        add_variables = true
        use_finite_deform_jacobian = true
        use_automatic_differentiation = true
      [../]
    [../]
  [../]
[]

[AuxVariables]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./TN_solid]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./TS_solid]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./T_x]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./T_y]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./T_z]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Tlocal_N]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Tlocal_S1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Tlocal_S2]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Functions]
  [./angles]
    type = PiecewiseLinear
    x = '0 1 2'
    y = '0 0 1.5707963'
  [../]

  [./stretch]
    type = PiecewiseLinear
    x = '0 1 2'
    y = '0 0.1 0.1'
  [../]

  [./move_y]
    type = ParsedFunction
    value = 'y*cos(theta) - z * (1 + a)*sin(theta) - y'
    vars = 'a theta'
    vals = 'stretch angles'
  [../]

  [./move_z]
    type = ParsedFunction
    value = 'y*sin(theta) + z*(1+a)*cos(theta) - z'
    vars = 'a theta'
    vals = 'stretch angles'
  [../]

  [./dt_fun]
    type = PiecewiseConstant
    x = '0 2'
    y = '0.05 0.05'
  []
[]

[BCs]
  [./fix]
    type = DirichletBC
    preset = true
    value = 0.0
    boundary = left
    variable = disp_x
  [../]

  [./front_y]
    type = FunctionDirichletBC
    boundary = front
    variable = disp_y
    function = move_y
    preset = true
  [../]

  [./back_y]
    type = FunctionDirichletBC
    boundary = back
    variable = disp_y
    function = move_y
    preset = true
  [../]

  [./front_z]
    type = FunctionDirichletBC
    boundary = front
    variable = disp_z
    function = move_z
    preset = true
  [../]

  [./back_z]
    type = FunctionDirichletBC
    boundary = back
    variable = disp_z
    function = move_z
    preset = true
  [../]

[]

[AuxKernels]
  # [./aux_Tavg_x]
  #   type = MaterialRealVectorValueAux
  #   boundary = 'interface'
  #   property = traction_global
  #   component = 0
  #   execute_on = 'TIMESTEP_END'
  #   variable = T_x
  # [../]
  # [./aux_Tavg_y]
  #   type = MaterialRealVectorValueAux
  #   boundary = 'interface'
  #   property = traction_global
  #   component = 1
  #   execute_on = 'TIMESTEP_END'
  #   variable = T_y
  # [../]
  # [./aux_Tavg_z]
  #   type = MaterialRealVectorValueAux
  #   boundary = 'interface'
  #   property = traction_global
  #   component = 2
  #   execute_on = 'TIMESTEP_END'
  #   variable = T_z
  # [../]
  # [./aux_Tloc_x]
  #   type = MaterialRealVectorValueAux
  #   boundary = 'interface'
  #   property = traction
  #   component = 0
  #   execute_on = 'TIMESTEP_END'
  #   variable = Tlocal_N
  # [../]
  # [./aux_Tloc_y]
  #   type = MaterialRealVectorValueAux
  #   boundary = 'interface'
  #   property = traction
  #   component = 1
  #   execute_on = 'TIMESTEP_END'
  #   variable = Tlocal_S1
  # [../]
  # [./aux_Tloc_z]
  #   type = MaterialRealVectorValueAux
  #   boundary = 'interface'
  #   property = traction
  #   component = 2
  #   execute_on = 'TIMESTEP_END'
  #   variable = Tlocal_S2
  # [../]
  [./stress_xx]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
    execute_on = timestep_end
    use_displaced_mesh = true
  [../]
  [./stress_yy]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
    execute_on = timestep_end
    use_displaced_mesh = true
  [../]
  [./stress_zz]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
    execute_on = timestep_end
    use_displaced_mesh = true
  [../]
  [./stress_xy]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
    execute_on = timestep_end
    use_displaced_mesh = true
  [../]
  [./stress_xz]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_xz
    index_i = 0
    index_j = 2
    execute_on = timestep_end
    use_displaced_mesh = true
  [../]
  [./stress_yz]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_yz
    index_i = 1
    index_j = 2
    execute_on = timestep_end
    use_displaced_mesh = true
  [../]
  # [./traction_N]
  #   type = TractionAux
  #   scalar_type = 'normal'
  #   variable = TN_solid
  #   property = stress
  #   boundary = 'interface'
  #   use_displaced_mesh = true
  # [../]
  # [./traction_S]
  #   type = TractionAux
  #   scalar_type = 'shear_norm'
  #   variable = TS_solid
  #   property = stress
  #   boundary = 'interface'
  #   use_displaced_mesh = true
  # [../]
[]



[Modules/TensorMechanics/CohesiveZoneMaster]
  [./czm]
    boundary = 'interface'
    displacements = 'disp_x disp_y disp_z'
  [../]
[]



# [Postprocessors]
#   [./sxx]
#     type = SideAverageValue
#     boundary = 'interface'
#     variable = stress_xx
#   [../]
#   [./syy]
#     type = SideAverageValue
#     boundary = 'interface'
#     variable = stress_yy
#   [../]
#   [./szz]
#     type = SideAverageValue
#     boundary = 'interface'
#     variable = stress_zz
#   [../]
#   [./syz]
#     type = SideAverageValue
#     boundary = 'interface'
#     variable = stress_yz
#   [../]
#   [./sxz]
#     type = SideAverageValue
#     boundary = 'interface'
#     variable = stress_xz
#   [../]
#   [./sxy]
#     type = SideAverageValue
#     boundary = 'interface'
#     variable = stress_xy
#   [../]
#   [./tn]
#     type = SideAverageValue
#     variable = TN_solid
#     boundary = 'interface'
#   [../]
#   [./ts]
#     type = SideAverageValue
#     variable = TS_solid
#     boundary = 'interface'
#   [../]
#   [./T_x]
#     type = SideAverageValue
#     variable = T_x
#     boundary = 'interface'
#   [../]
#   [./T_y]
#     type = SideAverageValue
#     variable = T_y
#     boundary = 'interface'
#   [../]
#   [./T_z]
#     type = SideAverageValue
#     variable = T_z
#     boundary = 'interface'
#   [../]
#   [./Tloc_N]
#     type = SideAverageValue
#     variable = Tlocal_N
#     boundary = 'interface'
#   [../]
#   [./Tloc_S1]
#     type = SideAverageValue
#     variable = Tlocal_S1
#     boundary = 'interface'
#   [../]
#   [./Tloc_S2]
#     type = SideAverageValue
#     variable = Tlocal_S2
#     boundary = 'interface'
#   [../]
# []

[Materials]
  [./stress]
    type = ADComputeFiniteStrainElasticStress
  [../]
  [./elasticity_tensor]
    type = ADComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '1.684e5 0.176e5 0.176e5 1.684e5 0.176e5 1.684e5 0.754e5 0.754e5 0.754e5'
  [../]
  [./czm]
    type = CZMMaterialLD
    displacements = 'disp_x disp_y disp_z'
    boundary = 'interface'
    large_kinematics = true
    E = 1e7
    G = 1e6
    use_area_change = true
    check_jacobian = false
  [../]
[]



[Preconditioning]
  [./smp]
    type = FDP
    full = true
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'newton'
  line_search = none

  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  l_max_its = 2
  l_tol = 1e-14
  nl_max_its = 15
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
  start_time = 0.0
  end_time = 2
  dtmin = 0.01
  [./TimeStepper]
    type = FunctionDT
    function = dt_fun
  [../]
[]

[Outputs]
  exodus = true
[]
