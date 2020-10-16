[Mesh]
  [./msh]
    type = FileMeshGenerator
    file = simple_diffusion_out_cp/0001_mesh.cpr
  []
  [./add_fake_neighbor]
    type = AssignFakeNeighborFromList
    input = msh
    fake_neighbor_list_file_name = 'simple_diffusion_bmbb.txt'
  []
[]

[Problem]
  #Note that the suffix is left off in the parameter below.
  restart_file_base = simple_diffusion_out_cp/0001  # You may also use a specific number here
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
  dt = 0.1
  end_time = 0.2
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
[]
