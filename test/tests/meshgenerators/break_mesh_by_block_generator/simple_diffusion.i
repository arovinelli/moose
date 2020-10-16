[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 2
  []
  [./subdomain_1]
    type = SubdomainBoundingBoxGenerator
    input = msh
    bottom_left = '0 0.5 0'
    top_right = '1 1 0'
    block_id = 1
  []
  [break]
    type = BreakMeshByBlockGenerator
    input = subdomain_1
    write_fake_neighbor_list_to_file = true
    fake_neighbor_list_file_name = 'simple_diffusion_bmbb.txt'
  []
[]

[Variables]
  [./u]
  [../]
[]

[Kernels]
  [./diff]
    type = CoefDiffusion
    variable = u
    coef = 0.1
  [../]
  [./time]
    type = TimeDerivative
    variable = u
  [../]
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  [../]
  [./right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 2
  dt = 0.1
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  [./out]
    type = Exodus
    sync_only = true
    sync_times = '0.1 0.2'
    execute_on = 'TIMESTEP_END'
  []
  checkpoint = true
[]
