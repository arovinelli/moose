#
# Stretch + rotation test
#
# This test is designed to compute a uniaxial stress and then follow it as the mesh is rotated by 45 degrees.
#
# The mesh is composed of two, single-elemnt blocks
# The large deforamtion traction separation kinematic assumes linear rotations and uses the velocity gradient L to keep track of area changes, hence it converges to the proper solutoin in the limit of dt->0. Smaller the time step higher the accuracy.

[Mesh]
  [./msh]
  type = GeneratedMeshGenerator
  dim = 3
  nx = 1
  ny = 1
  nz = 2
  xmin = -0.5
  xmax = 0.5
  ymin = -0.5
  ymax = 0.5
  zmin = -1
  zmax = 1
  []
  [./new_block]
    type = SubdomainBoundingBoxGenerator
    input = msh
    block_id = 1
    bottom_left = '-0.5 -0.5 0'
    top_right = '0.5 0.5 0.5'
  []
  [./split]
    type = BreakMeshByBlockGenerator
    input = new_block
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Functions]
  [./stretch]
    type = PiecewiseLinear
    x = '0 1'
    y = '0 1e2'
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
  [./back_z]
    type = FunctionNeumannBC
    boundary = front
    variable = disp_z
    use_displaced_mesh = false
    function = stretch
  [../]

[]


[Modules]
  [./TensorMechanics]
    [./Master]
      [./all]
        strain = FINITE
        add_variables = true
        use_finite_deform_jacobian = true
        use_automatic_differentiation = true
        generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_xz'
      [../]
    [../]
  [../]
[]

[Modules/TensorMechanics/CohesiveZoneMaster]
  [./czm1]
    boundary = 'interface'
    displacements = 'disp_x disp_y disp_z'
  [../]
[]

[Materials]
  [./stress]
    type = ADComputeFiniteStrainElasticStress
  [../]
  [./elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 1e3
    poissons_ratio = 0.3
  [../]
  [./czm_3dc]
    type = SalehaniIrani3DCTraction
    boundary = 'interface'
    normal_gap_at_maximum_normal_traction = 0.01
    tangential_gap_at_maximum_shear_traction = 0.005
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

[Postprocessors]
  [./Area]
     type = AreaPostprocessor
     boundary = interface
     use_displaced_mesh = true
  [../]
[]

[Executioner]
  # Executioner
  type = Transient

  solve_type = 'NEWTON'
  line_search = none
  petsc_options_iname = '-pc_type '
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-30
  nl_abs_tol = 1e-10
  l_max_its = 20
  start_time = 0.0
  dt = 0.25
  end_time = 0.25
[]

[Outputs]
  [./out]
    type = CSV
  [../]
[]
