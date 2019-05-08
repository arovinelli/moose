[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 2
  xmax = 2
  ny = 2
  ymax = 2
  nz = 2
  zmax = 2
[]

[MeshModifiers]
  [./subdomain1]
    type = SubdomainBoundingBox
    bottom_left = '0 0 0'
    top_right = '1 1 1'
    block_id = 1
  [../]
  [./break_boundary]
    depends_on = subdomain1
    type = BreakBoundaryOnSubdomain
  [../]
  [./interface]
    type = SideSetsBetweenSubdomains
    depends_on = break_boundary
    master_block = '0'
    paired_block = '1'
    new_boundary = 100
  [../]
  [./lower]
    depends_on = interface
    type = LowerDBlockFromSideset
    new_block_name = 'LD_interface'
    new_block_id = 100
    sidesets = '100'
  [../]
[]

[Variables]
  [./u]
    order = FIRST
    family = LAGRANGE
    block = 0
  [../]

  [./v]
    order = FIRST
    family = LAGRANGE
    block = 1
  [../]
  [./empty]
    block = 100
  [../]
[]

[Kernels]
  [./diff_u]
    type = CoeffParamDiffusion
    variable = u
    D = 4
    block = 0
  [../]
  [./diff_v]
    type = CoeffParamDiffusion
    variable = v
    D = 2
    block = 1
  [../]
  [./source_u]
    type = BodyForce
    variable = u
    value = 1
  [../]
  [./null]
    type = NullKernel
    variable = empty
  [../]
[]

[InterfaceKernels]
  [./interface]
    type = PenaltyInterfaceDiffusion
    variable = u
    neighbor_var = v
    boundary = 100
    penalty = 1e6
  [../]
[]

[UserObjects]
  [./map]
    type = map2LDelem
    boundary = 100
    execute_on = 'INITIAL'
    ld_block_names = 'LD_interface'
  [../]
[]

[BCs]
  [./u]
    type = VacuumBC
    variable = u
    boundary = 'left_to_0 bottom_to_0 back_to_0 right top front'
  [../]
  [./v]
    type = VacuumBC
    variable = v
    boundary = 'left_to_1 bottom_to_1 back_to_1'
  [../]
[]

[Postprocessors]
  [./u_int]
    type = ElementIntegralVariablePostprocessor
    variable = u
    block = 0
  [../]
  [./v_int]
    type = ElementIntegralVariablePostprocessor
    variable = v
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
  type = Steady
  solve_type = NEWTON
[]

[Outputs]
  exodus = true
  print_linear_residuals = true
[]
