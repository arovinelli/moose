[Mesh]
  [msh]
    type = FileMeshGenerator
    file = split_FMG_diffusion.cpr
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
