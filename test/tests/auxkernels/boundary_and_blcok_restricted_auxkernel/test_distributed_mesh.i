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
      subdomain_ids = '0 1
                       1 1'
    [../]

    [./interface]
      type = SideSetsBetweenSubdomains
      depends_on = subdomain_id
      master_block = '0'
      paired_block = '1'
      new_boundary = 'interface'
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
  []

  [AuxKernels]
    [./diffusivity_1]
      type = MaterialRealAux
      property = diffusivity
      variable = diffusivity_1
      block = 0
    []
    [./diffusivity_2]
      type = MaterialRealAux
      property = diffusivity
      variable = diffusivity_2
      block = 1
    []
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
  []


  [Executioner]
    type = Steady
    solve_type = NEWTON
  []

  [Outputs]
    exodus = true
  []
