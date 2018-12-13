[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = -1
  ymin = -1
  xmax = 1
  ymax = 1
  nx = 2
  ny = 2
  elem_type = QUAD4
[]

[MeshModifiers]
  [./subdomain_id]
    type = AssignElementSubdomainID
    subdomain_ids = '0 0
                     1 2'
  [../]

  [./interface_01]
    type = SideSetsBetweenSubdomains
    depends_on = subdomain_id
    master_block = '0'
    paired_block = '1'
    new_boundary = 'interface_01'
  [../]
  [./interface_02]
    type = SideSetsBetweenSubdomains
    depends_on = subdomain_id
    master_block = '0'
    paired_block = '2'
    new_boundary = 'interface_02'
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
  [./interface_uo_qp]
    type = InterfaceUO_QP
    variable = u
    diffusivity = diffusivity
    boundary = 'interface_01 interface_02'
    execute_on = 'INITIAL LINEAR'
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
  [./stateful1]
    type = StatefulMaterial
    block = 0
    initial_diffusivity = 5
  [../]
  [./stateful2]
    type = StatefulMaterial
    block = 1
    initial_diffusivity = 2
  [../]
  [./stateful3]
    type = StatefulMaterial
    block = 2
    initial_diffusivity = 7
  [../]
  [./interface_material]
    type = InterfaceUOMaterial
    is_interface_material = true
    boundary = 'interface_01 interface_02'
    interface_uo_qp = interface_uo_qp
  [../]
[]

[AuxVariables]
  [./ujump_01]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./ujump_02]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./boundary_property_01]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./boundary_property_02]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]


[AuxKernels]
  [./ujump_01]
    type = MaterialRealAux
    property = variable_jump
    variable = ujump_01
    boundary = 'interface_01'
  [../]
  [./ujump_02]
    type = MaterialRealAux
    property = variable_jump
    variable = ujump_02
    boundary = 'interface_02'
  [../]
  [./boundary_property_01]
    type = MaterialRealAux
    property = boundary_property
    variable = boundary_property_01
    boundary = 'interface_01'
  [../]
  [./boundary_property_02]
    type = MaterialRealAux
    property = boundary_property
    variable = boundary_property_02
    boundary = 'interface_02'
  [../]
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
[]

[Outputs]
  exodus = true
[]
