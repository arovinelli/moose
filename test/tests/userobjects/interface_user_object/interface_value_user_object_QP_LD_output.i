[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 2
  xmax = 2
  ny = 2
  ymax = 2
  elem_type = QUAD4
[]

[MeshModifiers]
  [./subdomain_id1]
    type = SubdomainBoundingBox
    bottom_left = '0 0 0'
    top_right = '1 1 0'
    block_id = 1
  [../]
  [./subdomain_id2]
    depends_on = subdomain_id1
    type = SubdomainBoundingBox
    bottom_left = '1 0 0'
    top_right = '2 1 0'
    block_id = 2
  [../]

  [./interface]
    type = SideSetsBetweenSubdomains
    depends_on = subdomain_id2
    master_block = '0 0 1'
    paired_block = '1 2 2'
    new_boundary = '100'
  [../]

  [./lower]
    depends_on = interface
    type = LowerDBlockFromSideset
    new_block_name = 'LD_interface'
    new_block_id = 50
    sidesets = '100'
  [../]

[]

[Functions]
  [./fn_exact]
    type = ParsedFunction
    value = 'x*x+y*y'
  [../]

  [./ffn]
    type = ParsedFunction
    value = -4
  [../]
[]

[UserObjects]
  [./interface_value_uo]
    type = map2LDelem
    var = diffusivity_1
    var_neighbor = diffusivity_2
    boundary = '100'
    execute_on = 'initial timestep_end'
    interface_value_type = average
    ld_block_names = 'LD_interface'
  [../]
  [./interface_master_minus_slave_uo]
    type = map2LDelem
    var = diffusivity_1
    var_neighbor = diffusivity_2
    boundary = '100'
    execute_on = 'initial timestep_end'
    interface_value_type = jump_master_minus_slave
    ld_block_names = 'LD_interface'
  [../]
  [./interface_slave_minus_master_uo]
    type = map2LDelem
    var = diffusivity_1
    var_neighbor = diffusivity_2
    boundary = '100'
    execute_on = 'initial timestep_end'
    interface_value_type = jump_slave_minus_master
    ld_block_names = 'LD_interface'
  [../]
  [./interface_absolute_jump_uo]
    type = map2LDelem
    var = diffusivity_1
    var_neighbor = diffusivity_2
    boundary = '100'
    execute_on = 'initial timestep_end'
    interface_value_type = jump_abs
    ld_block_names = 'LD_interface'
  [../]
  [./interface_master_uo]
    type = map2LDelem
    var = diffusivity_1
    var_neighbor = diffusivity_2
    boundary = '100'
    execute_on = 'initial timestep_end'
    interface_value_type = master
    ld_block_names = 'LD_interface'
  [../]
  [./interface_slave_uo]
    type = map2LDelem
    var = diffusivity_1
    var_neighbor = diffusivity_2
    boundary = '100'
    execute_on = 'initial timestep_end'
    interface_value_type = slave
    ld_block_names = 'LD_interface'
  [../]
[]

[Variables]
  [./u]
    family = LAGRANGE
    order = FIRST
  [../]
[]


[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]

  [./ffn]
    type = BodyForce
    variable = u
    function = ffn
  [../]
[]

[BCs]
  [./all]
    type = FunctionDirichletBC
    variable = u
    boundary = '0 1 2 3'
    function = fn_exact
  [../]
[]

[Materials]
  [./stateful0]
    type = StatefulMaterial
    block = 0
    initial_diffusivity = 3
    # outputs = all
  [../]
  [./stateful1]
    type = StatefulMaterial
    block = 1
    initial_diffusivity = 6
    # outputs = all
  [../]
  [./stateful2]
    type = StatefulMaterial
    block = 2
    initial_diffusivity = 12
    # outputs = all
  [../]
  [./dummy]
    type = GenericConstantMaterial
    block = 50
  [../]
[]

[AuxKernels]
  [./diffusivity_1]
    type = MaterialRealAux
    property = diffusivity
    variable = diffusivity_1
    block = '0 1 2'
  []
  [./diffusivity_2]
    type = MaterialRealAux
    property = diffusivity
    variable = diffusivity_2
    block = '0 1 2'
  []
  [./interface_avg_qp_aux]
    type = InterfaceValueUserObjectAuxLD
    variable = avg_qp
    block = '50'
    interface_uo_name = interface_value_uo
  []
  [./interface_master_minus_slave_qp_aux]
    type = InterfaceValueUserObjectAuxLD
    variable = master_minus_slave_qp
    block = '50'
    interface_uo_name = interface_master_minus_slave_uo
  [../]
  [./interface_slave_minus_master_qp_aux]
    type = InterfaceValueUserObjectAuxLD
    variable = slave_minus_master_qp
    block = '50'
    interface_uo_name = interface_slave_minus_master_uo
  [../]
  [./interface_absolute_jump_qp_aux]
    type = InterfaceValueUserObjectAuxLD
    variable = abs_jump_qp
    block = '50'
    interface_uo_name = interface_absolute_jump_uo
  [../]
  [./interface_master_qp_aux]
    type = InterfaceValueUserObjectAuxLD
    variable = master_qp
    block = '50'
    interface_uo_name = interface_master_uo
  [../]
  [./interface_slave_qp_aux]
    type = InterfaceValueUserObjectAuxLD
    variable = slave_qp
    block = '50'
    interface_uo_name = interface_slave_uo
  [../]
[]

[AuxVariables]
  [./diffusivity_1]
    family = MONOMIAL
    order = CONSTANT
  []
  [./diffusivity_2]
    family = MONOMIAL
    order = CONSTANT
  []
  [./avg_qp]
    family = MONOMIAL
    order = CONSTANT
  []
  [./master_minus_slave_qp]
    family = MONOMIAL
    order = CONSTANT
  []
  [./slave_minus_master_qp]
    family = MONOMIAL
    order = CONSTANT
  []
  [./abs_jump_qp]
    family = MONOMIAL
    order = CONSTANT
  []
  [./master_qp]
    family = MONOMIAL
    order = CONSTANT
  []
  [./slave_qp]
    family = MONOMIAL
    order = CONSTANT
  []
[]



[Postprocessors]
  [./interface_average_PP]
    type = ElementAverageValue
    block = 50
    variable =  avg_qp
  [../]
  [./master_minus_slave_qp_PP]
    type = ElementAverageValue
    block = 50
    variable =  master_minus_slave_qp
  [../]
  [./slave_minus_master_qp_PP]
    type = ElementAverageValue
    block = 50
    variable =  slave_minus_master_qp
  [../]
  [./abs_jump_qp_PP]
    type = ElementAverageValue
    block = 50
    variable =  abs_jump_qp
  [../]
  [./master_qp_PP]
    type = ElementAverageValue
    block = 50
    variable =  master_qp
  [../]
  [./slave_qp_PP]
    type = ElementAverageValue
    block = 50
    variable =  slave_qp
  [../]
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
[]

[Outputs]
  exodus = true
[]
