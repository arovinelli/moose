
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
      [./u]
        [./InitialCondition]
          type = ConstantIC
          value = 1
        [../]
        block = 0
      [../]
      [./u_neighbor]
        [./InitialCondition]
          type = ConstantIC
          value = 1
        [../]
        block = 1
      [../]
[]

[Mesh]
  [./msh]
  type = GeneratedMeshGenerator
  dim = 3
  nx = 1
  ny = 1
  nz = 2
  xmin = -1
  xmax = 0
  zmin = -1
  zmax = 1
  []
  [./new_block]
    type = SubdomainBoundingBoxGenerator
    input = msh
    block_id = 1
    bottom_left = '-1 0 0'
    top_right = '0 1 1'
  []
  [./split]
    type = BreakMeshByBlockGenerator
    input = new_block
  []
  [./top_nodes_left]
    input = split
    type = ExtraNodesetGenerator
    nodes = '8 11 12 15'
    new_boundary = top_nodes_left
  []
  [./top_nodes_right]
    input = top_nodes_left
    type = ExtraNodesetGenerator
    nodes = '9 10 13 14'
    new_boundary = top_nodes_right
  []
  [./bottom_nodes_left]
    input = top_nodes_right
    type = ExtraNodesetGenerator
    nodes = '0 3 4 7'
    new_boundary = bottom_nodes_left
  []
  [./bottom_nodes_right]
    input = bottom_nodes_left
    type = ExtraNodesetGenerator
    nodes = '1 2 5 6 '
    new_boundary = bottom_nodes_right
  []
[]

[InterfaceKernels]
  [./flux_match]
    type = PenaltyInterfaceDiffusion
    variable = u
    neighbor_var = u_neighbor
    boundary = interface
    penalty = 1e6
    use_displaced_mesh = true
  [../]
[]

[Materials]
  [./stateful]
    type = StatefulTest
    prop_names = 'diffusivity'
    prop_values = '1'
    block = '0 1'
  [../]
[]

[Kernels]
  [nullx]
    type =  NullKernel
    variable = disp_x
  []
  [nully]
    type =  NullKernel
    variable = disp_y
  []
  [nullz]
    type =  NullKernel
    variable = disp_z
  []
  [./diff]
    type = MatDiffusionTest
    variable = u
    prop_name = diffusivity
    block = 0
  [../]
  [./abs]
    type = Reaction
    variable = u
    block = 0
  [../]
  [./forcing]
    type = BodyForce
    variable = u
    function = forcing_fn
    block = 0
  [../]
  [./diffn]
    type = MatDiffusionTest
    variable = u_neighbor
    prop_name = diffusivity
    block = 1
  [../]
  [./absn]
    type = Reaction
    variable = u_neighbor
    block = 1
  [../]
  [./forcingn]
    type = BodyForce
    variable = u_neighbor
    function = forcing_fn
    block = 1
  [../]
[]


[Preconditioning]
  [./smp]
    type = SMP
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
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6
  dtmin = 1
  dtmax = 1
  end_time = 2.5
[]

[Functions]

  [./move_x]
    type = ParsedFunction
    value = 0.1*t
  [../]
  [./move_x2]
    type = ParsedFunction
    value = -0.1*t
  [../]
[]


[BCs]
  [./u]
    type = FunctionDirichletBC
    variable = u
    boundary = 'left'
    function = bc_fn
  [../]
  [./u_neighbor]
    type = FunctionDirichletBC
    variable = u_neighbor
    boundary = 'right'
    function = bc_fn
  [../]
  [./x_move]
    type = FunctionDirichletBC
    boundary = 'top_nodes_right'
    variable = disp_x
    function = move_x
  [../]
  [./x_move_2]
    type = FunctionDirichletBC
    boundary = 'bottom_nodes_right'
    variable = disp_x
    function = move_x2
  [../]


  [./x_fix]
    type = DirichletBC
    boundary = 'bottom_nodes_left top_nodes_left'
    variable = disp_x
    value = 0.0
  [../]
  [./y_fix]
    type = DirichletBC
    boundary = 'bottom_nodes_left bottom_nodes_right top_nodes_left top_nodes_right '
    variable = disp_y
    value = 0.0
  [../]
  [./z_fix]
    type = DirichletBC
    boundary = 'bottom_nodes_left bottom_nodes_right top_nodes_left top_nodes_right'
    variable = disp_z
    value = 0.0
  [../]
[]

[Functions]
  [./forcing_fn]
    type = ParsedFunction
    value = (x*x*x)-6.0*x
  [../]

  [./bc_fn]
    type = ParsedFunction
    value = (x*x*x)
  [../]
[]


[Outputs]
  exodus = true
  sync_times = '0 0.5 1 1.5 2'
[]
