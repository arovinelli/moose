[Mesh]
  file = czm_non_local_in.e
  parallel_type = REPLICATED
[]

[MeshModifiers]
  # [./subdomain1]
  #   type = SubdomainBoundingBox
  #   bottom_left = '0 0 0'
  #   top_right = '1 1 1'
  #   block_id = 1
  # [../]
  # [./subdomain2]
  #   type = SubdomainBoundingBox
  #   bottom_left = '0 0 1'
  #   top_right = '1 1 2'
  #   block_id = 2
  #   depends_on = subdomain1
  # [../]

  [./breakmesh]
    type = BreakMeshByBlock
  [../]
  #
  # [./bottom_block_1]
  #   type = SideSetsAroundSubdomain
  #   depends_on = 'breakmesh'
  #   block = '1'
  #   new_boundary = 'back'
  #   normal = '0 0 -1'
  # [../]
  # [./top_block_2]
  #   type = SideSetsAroundSubdomain
  #   depends_on = 'breakmesh'
  #   block = '2'
  #   new_boundary = 'front'
  #   normal = '0 0 1'
  # [../]
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]



[Modules/TensorMechanics/Master]
  [./all]
    strain = SMALL
    add_variables = true
    generate_output = 'strain_yy strain_zz strain_yz strain_xz strain_xy'
  [../]
[]


[BCs]
  [./bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = back
    value = 0.0
  [../]
  [./bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = back
    value = 0.0
  [../]
  [./bottom_z]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]
  [./top2_x]
    type = DirichletBC
    variable = disp_x
    boundary = front
    value = 0.0
  [../]
  [./top2_y]
    type = DirichletBC
    variable = disp_y
    boundary = front
    value = 0.0
  [../]
  [./top2_z]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = front
    function = 1*t
  [../]
[]
[InterfaceKernels]
  [./interface_x]
    # type = CZMInterfaceKernel
    type = CZMInterfaceKernelWithAvgAuxVars
    variable = disp_x
    neighbor_var = disp_x
    disp_1 = disp_y
    disp_1_neighbor = disp_y
    disp_2 = disp_z
    disp_2_neighbor = disp_z
    disp_index = 0
    boundary = 'interface'
    var_0_name = avg_stress_vm_interface
  [../]
  [./interface_y]
    # type = CZMInterfaceKernel
    type = CZMInterfaceKernelWithAvgAuxVars
    variable = disp_y
    neighbor_var = disp_y
    disp_1 = disp_x
    disp_1_neighbor = disp_x
    disp_2 = disp_z
    disp_2_neighbor = disp_z
    disp_index = 1
    boundary = 'interface'
    var_0_name = avg_stress_vm_interface
  [../]
  [./interface_z]
    # type = CZMInterfaceKernel
    type = CZMInterfaceKernelWithAvgAuxVars
    variable = disp_z
    neighbor_var = disp_z
    disp_1 = disp_x
    disp_1_neighbor = disp_x
    disp_2 = disp_y
    disp_2_neighbor = disp_y
    disp_index = 2
    boundary = 'interface'
    var_0_name = avg_stress_vm_interface
  [../]
[]
[AuxVariables]
  [./vonmises_stress]
    family = MONOMIAL
    order = CONSTANT
    [./InitialCondition]
      type = ConstantIC
      value = 0.
    [../]
  []
  # [./stress_zz]
  #   family = MONOMIAL
  #   order = CONSTANT
  # []
  # [./stress_hyd]
  #   family = MONOMIAL
  #   order = CONSTANT
  # []
  # [./avg_stress_zz_interface]
  #   family = MONOMIAL
  #   order = CONSTANT
  # []
  [./avg_stress_vm_interface]
    family = MONOMIAL
    order = CONSTANT
  []
  [./avg_jump_z]
    family = MONOMIAL
    order = CONSTANT
  []
[]


[AuxKernels]
  [./stress_vm_aux]
    type = RankTwoScalarAux
    scalar_type = VonMisesStress
    rank_two_tensor = stress
    variable = vonmises_stress
    execute_on = 'timestep_end'
  []
  # [./stress_zz_aux]
  #   type = RankTwoAux
  #   index_i = 2
  #   index_j = 2
  #   rank_two_tensor = stress
  #   variable = stress_zz
  #   execute_on = 'initial NONLINEAR LINEAR timestep_end timestep_begin'
  # []
  [./jump_z_aux]
    type = MaterialRealVectorValueAux
    component = 2
    variable = avg_jump_z
    property = displacement_jump
    execute_on = ' NONLINEAR LINEAR timestep_end timestep_begin'
    boundary = interface
  []
  # [./stress_hyd_aux]
  #   type = RankTwoScalarAux
  #   scalar_type = Hydrostatic
  #   rank_two_tensor = stress
  #   variable = stress_hyd
  # []
  # [./stress_zz_interface]
  #   type = MaterialRealAux
  #   property = stress_zz_mp
  #   variable = avg_stress_zz_interface
  #   boundary = interface
  #   execute_on = 'NONLINEAR LINEAR timestep_end timestep_begin'
  # []
  [./stress_vm_interface]
    type = MaterialRealAux
    property = stress_vm_mp
    variable = avg_stress_vm_interface
    boundary = interface
    execute_on = 'NONLINEAR LINEAR timestep_end'
  []
  # [./stress_hyd_interface]
  #   type = MaterialRealAux
  #   property = stress_hyd_mp
  #   variable = avg_stress_hyd_interface
  #   boundary = interface
  # []
[]

[UserObjects]
  # [./avg_interface_stress_zz]
  #   type = ScalarBulkMPAcrossInterface_QP
  #   var_name = stress_zz
  #   boundary = 'interface'
  #   execute_on = 'initial NONLINEAR LINEAR timestep_end timestep_begin'
  # []
  [./avg_interface_stress_vm]
    type = ScalarBulkMPAcrossInterface_QP
    var_name = vonmises_stress
    boundary = 'interface'
    execute_on = 'initial timestep_begin'
    # use_old_prop = true
  []
  # [./avg_interface_hydstress]
  #   type = ScalarBulkMPAcrossInterface_QP
  #   var_name = hydrostatic_stress
  #   boundary = 'interface'
  #   execute_on = 'initial NONLINEAR LINEAR timestep_end'
  # []
  [./displacement_jump_uo]
    type = DispJumpUO_QP
    disp_x = disp_x
    disp_y = disp_y
    disp_z = disp_z
    boundary = 'interface'
    execute_on = 'initial NONLINEAR LINEAR timestep_end timestep_begin'
  [../]
  [./cohesive_law_SVM]
    type = CZMTestLaw_SVM
    # MaxAllowableTraction = '100 70'
    # DeltaU0 = '1 0.7'
    displacement_jump_mp_name = 'displacement_jump_local'
    boundary = 'interface'
    other_avg_scalar_mp_names = 'stress_vm_mp'
    n_other_avg_scalar_mp_names = 1
    interface_stiffness = 0.5e0
  [../]
[]

[Materials]
  [./Elasticity_tensor_b1]
    type = ComputeElasticityTensor
    block = '1'
    fill_method = symmetric_isotropic
    C_ijkl = '0.3 0.5e8'
  [../]
  [./Elasticity_tensor_b2]
    type = ComputeElasticityTensor
    block = '2'
    fill_method = symmetric_isotropic
    C_ijkl = '0.3 0.5e7'
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = '1 2'
  [../]
  [./gap]
    type = CZMUOBasedMaterial
    is_interface_material = true
    boundary = 'interface'
    displacement_jump_UO = 'displacement_jump_uo'
    traction_separation_UO = 'cohesive_law_SVM'
    bulk_avg_mp_names = 'stress_vm_mp'
    bulk_avg_mp_uo_names = 'avg_interface_stress_vm'
    need_avg_scalar_vars_derivatives = true
  [../]
[]
 [Preconditioning]
   [./SMP]
     type = SMP
     full = true
   [../]
 []
[Executioner]
  # Preconditisoned JFNK (default)
  type = Transient
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  # petsc_options_value = 'hypre     boomerang'
  solve_type = NEWTON
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-6
  nl_max_its = 20
  l_tol = 1e-10
  l_max_its = 50
  start_time = 0.0
  dt = 5
  end_time = 100
  # dtmin = 1
  line_search = none
  # num_steps = 1
[]
[Outputs]
  [./out]
    type = Exodus
  [../]
[]
[Postprocessors]
  # [./sxx_3G]
  #   type = ElementAverageValue
  #   variable = stress_xx
  #   execute_on = 'initial timestep_end'
  #   block = 2
  # [../]
  # [./syy_3G]
  #   type = ElementAverageValue
  #   variable = stress_yy
  #   execute_on = 'initial timestep_end'
  #   block = 2
  # [../]
  # [./szz_3G]
  #   type = ElementAverageValue
  #   variable = stress_zz
  #   execute_on = 'initial timestep_end'
  #   block = 2
  # [../]
  [./svm_3G]
    type = ElementAverageValue
    variable = vonmises_stress
    execute_on = 'initial timestep_end'
    block = 2
  [../]
  # [./syz_3G]
  #   type = ElementAverageValue
  #   variable = stress_yz
  #   execute_on = 'initial timestep_end'
  #   block = 2
  # [../]
  # [./sxz_3G]
  #   type = ElementAverageValue
  #   variable = stress_xz
  #   execute_on = 'initial timestep_end'
  #   block = 2
  # [../]
  # [./sxy_3G]
  #   type = ElementAverageValue
  #   variable = stress_xy
  #   execute_on = 'initial timestep_end'
  #   block = 2
  # [../]
  [./disp_top3_z]
    type = SideAverageValue
    variable = disp_z
    execute_on = 'initial timestep_end'
    boundary = 'front'
  [../]
[]
