# Patch test for cohesive zone modeling to check convergence
[Mesh]
  [./msh]
    type = FileMeshGenerator
    file = czm_patch_test_only_stretch_reference_out_cp/0003_mesh.cpr
  []
  [./add_fake_neighbor]
    type = AssignFakeNeighborFromList
    input = msh
  []
[]

[Problem]
  #Note that the suffix is left off in the parameter below.
  restart_file_base = czm_patch_test_only_stretch_reference_out_cp/0003  # You may also use a specific number here
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Modules]
  [./TensorMechanics]
    [./Master]
      [./all]
        strain = FINITE
        add_variables = true
        use_finite_deform_jacobian = true
        use_automatic_differentiation = true
        generate_output = 'stress_xx stress_yy stress_zz stress_yz stress_xz stress_xy'
      [../]
    [../]
  [../]
[]


[Functions]
  [./stretch]
    type = PiecewiseLinear
    x = '0 1'
    y = '0 0.01'
  [../]
[]

[BCs]
  [./fix_x]
    type = DirichletBC
    preset = true
    value = 0.0
    boundary = left
    variable = disp_x
  [../]
  [./fix_y]
    type = DirichletBC
    preset = true
    value = 0.0
    boundary = bottom
    variable = disp_y
  [../]
  [./fix_z]
    type = DirichletBC
    preset = true
    value = 0.0
    boundary = back
    variable = disp_z
  [../]


  [./stretch_z]
    type = FunctionDirichletBC
    boundary = front
    variable = disp_z
    function = stretch
    preset = true
  [../]
[]


[Modules/TensorMechanics/CohesiveZoneMaster]
  [./czm_ik]
    boundary = 'interface'
  [../]
[]

[Materials]
  [./stress]
    type = ADComputeFiniteStrainElasticStress
  [../]
  [./elasticity_tensor]
    type = ADComputeElasticityTensor
    fill_method = symmetric9
    C_ijkl = '1.684e5 0.176e5 0.176e5 1.684e5 0.176e5 1.684e5 0.754e5 0.754e5 0.754e5'
  [../]
  [./czm_mat]
    boundary = 'interface'
    type=SalehaniIrani3DCTraction
    normal_gap_at_maximum_normal_traction=0.05
    tangential_gap_at_maximum_shear_traction=0.05
    maximum_normal_traction=1e1
    maximum_shear_traction=7e0
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
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
  dt = 0.25
  end_time = 1
[]

[Postprocessors]
  [./nonlin]
    type = NumNonlinearIterations
  [../]
[]


[Outputs]
  checkpoint = true
  [./out]
    type = Exodus
    sync_only = true
    sync_times = '1'
    execute_on = 'TIMESTEP_END'
  []
[]
