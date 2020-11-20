[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    nx=10
    ny=10
    nz=10
    elem_type = HEX8
    dim=3
  []
  [add_side_sets]
    input = msh
    type = SideSetsFromNormalsGenerator
    normals = '0 -1  0
               0  1  0
               -1 0  0
               1  0  0
               0  0 -1
               0  0  1'
    fixed_normal = true
    new_boundary = 'y0 y1 x0 x1 z0 z1'
  []
[]


[Variables]
  [u]
  []
[]

[Functions]
  [top_bc]
    type = ParsedFunction
    value = 'x'
  []
[]

[Kernels]
  [diff]
    type = AnisotropicDiffusion
    variable = u
    tensor_coeff = '2 0 0
                    0 4 0
        0 0 6'
  []
  [time]
    type = TimeDerivative
    variable = u
  []
[]

[BCs]
  [lower_left]
    type = DirichletBC
    variable = u
    boundary = 'y0 x0'
    value = 1
  []

  [top]
    type = FunctionNeumannBC
    variable = u
    boundary = y1
    function = top_bc
  []

  [right]
    type = NeumannBC
    variable = u
    boundary = x1
  []
[]

[Executioner]
  type = Transient
  num_steps = 2
  dt = 0.1
  solve_type = PJFNK
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]
