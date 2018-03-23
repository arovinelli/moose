[Mesh]
  file = testInputTractionNormals.e
    displacements = 'disp_x disp_y disp_z'
  [./MortarInterfaces]
    [./middle]
      master = 100
      slave = 101
      subdomain = 1000
    [../]
  [../]
[]

# [GlobalParams]
#      use_displaced_mesh = true
# []

[NodalNormals]
  boundary = '100 101'
  corner_boundary = '202'
  order = FIRST
[]

# [NodalNormals]
#   boundary = '101'
#   corner_boundary = 301
#   order = SECOND
# []


[Problem]
      coord_type = XYZ
[]

[Functions]
  [./rampDispl]
     type = PiecewiseLinear
     x = '0 1'
     y = '0 0.1'
  [../]

  [./jumpX]
      type = ParsedFunction
      value = t*0.01
  [../]

  [./jumpY]
      type = ParsedFunction
      value = t*0.02
  [../]

  [./jumpZ]
      type = ParsedFunction
      value = t*0.5
  [../]


  [./damaget]
      type = ParsedFunction
      value = 1./(1.-t)*1e-8
  [../]

  [./jump0]
      type = ParsedFunction
      value = t*0.
  [../]

[]

[Variables]
      [./disp_x]
        block = '1 2'
      [../]
      [./disp_y]
        block = '1 2'
      [../]
      [./disp_z]
        block = '1 2'
      [../]
  [./ux_lm]
    order = FIRST
    family = LAGRANGE
    block = middle
  [../]
  [./uy_lm]
    order = FIRST
    family = LAGRANGE
    block = middle
  [../]
  [./uz_lm]
    order = FIRST
    family = LAGRANGE
    block = middle
  [../]
[]

[AuxVariables]
  [./tX]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tY]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tZ]
    order = CONSTANT
    family = MONOMIAL
  [../]

  # [./Ntraction]
  #   order = FIRST
  #   family = MONOMIAL
  # [../]

  [./tFX]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tFY]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tFZ]
    order = CONSTANT
    family = MONOMIAL
  [../]


  [./s00]
    order = FIRST
    family = MONOMIAL
  [../]
  [./s01]
    order = FIRST
    family = MONOMIAL
  [../]
  [./s02]
    order = FIRST
    family = MONOMIAL
  [../]
  [./s10]
    order = FIRST
    family = MONOMIAL
  [../]
  [./s11]
    order = FIRST
    family = MONOMIAL
  [../]
  [./s12]
    order = FIRST
    family = MONOMIAL
  [../]
  [./s20]
    order = FIRST
    family = MONOMIAL
  [../]
  [./s21]
    order = FIRST
    family = MONOMIAL
  [../]
  [./s22]
    order = FIRST
    family = MONOMIAL
  [../]

  [./NodalForce_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./NodalForce_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./NodalForce_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxKernels]
  [./s00]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    variable = s00
    # boundary = '100 101'
    block = '1 2 '
  []
  [./s01]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    variable = s01
    # boundary = '100 101'
    block = '1 2 '
  []
  [./s02]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 2
    variable = s02
    # boundary = '100 101'
    block = '1 2 '
  []
  [./s10]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 0
    variable = s10
    # boundary = '100 101'
    block = '1 2 '
  []
  [./s11]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    variable = s11
    # boundary = '100 101'
    block = '1 2 '
  []
  [./s12]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 2
    variable = s12
    # boundary = '100 101'
    block = '1 2 '
  []
  [./s20]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 2
    index_j = 0
    variable = s20
    # boundary = '100 101'
    block = '1 2 '
  []
  [./s21]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 2
    index_j = 1
    variable = s21
    # boundary = '100 101'
    block = '1 2 '
  []
  [./s22]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 2
    index_j = 2
    variable = s22
    # boundary = '100 101'
    block = '1 2 '
  []
  [./trationX]
    type = TractionAux
    stress= stress
    tIdx = 0
    variable = tX
    boundary = '100 101'
  [../]
  [./trationY]
    type = TractionAux
    stress= stress
    tIdx = 1
    variable = tY
    boundary = '100 101'
  [../]
  [./trationZ]
    type = TractionAux
    stress= stress
    tIdx = 2
    variable = tZ
    boundary = '100 101'
  [../]
  # [./tractionN]
  #   type = NormalTractionAux
  #   stress= stress
  #   boundary = '100 101'
  #   variable =  Ntraction
  # [../]

  [./trationXF]
    type = TractionAuxField
    s00= s00
    s01= s01
    s02= s02
    s10= s10
    s11= s11
    s12= s12
    s20= s20
    s21= s21
    s22= s22
    tIdx = 0
    variable = tFX
    boundary = '100 101'
  [../]
  [./trationYF]
    type = TractionAuxField
    s00= s00
    s01= s01
    s02= s02
    s10= s10
    s11= s11
    s12= s12
    s20= s20
    s21= s21
    s22= s22
    tIdx = 1
    variable = tFY
    boundary = '100 101'
  [../]
  [./trationZF]
    type = TractionAuxField
    s00= s00
    s01= s01
    s02= s02
    s10= s10
    s11= s11
    s12= s12
    s20= s20
    s21= s21
    s22= s22
    tIdx = 2
    variable = tFZ
    boundary = '100 101'
  [../]
 []


[Kernels]
  [./TensorMechanics]
     use_displaced_mesh = true
     block= '1 2'
     displacements = 'disp_x disp_y disp_z'
     save_in = 'NodalForce_x NodalForce_y NodalForce_z'
  [../]
[]

[Constraints]
  [./UX]
    type = EqualValueConstraintWithJumpFuncANL
    variable = ux_lm
    interface = middle
    # traction = tX
    funcJump= jumpX
    # damageFunc= damaget
    master_variable = disp_x
  [../]
  [./UY]
    type = EqualValueConstraintWithJumpFuncANL
    variable = uy_lm
    interface = middle
    funcJump= jump0
    master_variable = disp_y
  [../]
  [./UZ]
    type = EqualValueConstraintWithJumpFuncANL
    variable = uz_lm
    interface = middle
    funcJump= jump0
    master_variable = disp_z
  [../]
[]



[Materials]
  [./elasticity_tensor1]
    type = ComputeIsotropicElasticityTensor
    block = '1'
    youngs_modulus = 100e3
    poissons_ratio= 0.3
  [../]

  [./elasticity_tensor2]
    type = ComputeIsotropicElasticityTensor
    block = '2'
    youngs_modulus = 100e3
    poissons_ratio= 0.3
  [../]

  [./strain]
    type = ComputeFiniteStrain
    displacements = 'disp_x disp_y disp_z'
    block = '1 2'
  [../]

  [./finite_strain_elastic_stress]
    type = ComputeFiniteStrainElasticStress
    block = '1 2'
  [../]
[]



[BCs]



      [./nox]
            type = PresetBC
            variable = disp_x
            boundary = '1'
            value = 0.0
      [../]
      [./noy]
            type = PresetBC
            variable = disp_y
            boundary = '1'
            value = 0.0
      [../]
      [./noz]
            type = PresetBC
            variable = disp_z
            boundary = '1'
            value = 0.0
      [../]



      [./Loads]
            type = FunctionPresetBC
            variable = disp_x
            boundary = '2'
            function = rampDispl
      [../]


[]


[Preconditioning]
  [./fmp]
    type = SMP
    full = true

  [../]
[]

[Executioner]
    solve_type = 'NEWTON'
  type = Transient
  l_max_its = 500
#  l_tol = 1e-8
#  nl_max_its = 35
#  nl_rel_tol = 1e-6
#  nl_abs_tol = 1e-6

#  petsc_options = '-snes_converged_reason -ksp_converged_reason'
#  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
#  petsc_options_value = 'lu superlu_dist'
  dt = 0.1
  end_time = 1
[]

[Outputs]

  exodus = true
[]
