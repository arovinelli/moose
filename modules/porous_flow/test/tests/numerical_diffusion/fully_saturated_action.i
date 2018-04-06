# Using the fully-saturated action, which does mass lumping but no upwinding
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 100
  xmin = 0
  xmax = 1
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [./porepressure]
  [../]
  [./tracer]
  [../]
[]

[ICs]
  [./porepressure]
    type = FunctionIC
    variable = porepressure
    function = '1 - x'
  [../]
  [./tracer]
    type = FunctionIC
    variable = tracer
    function = 'if(x<0.1,0,if(x>0.3,0,1))'
  [../]
[]

[PorousFlowFullySaturated]
  porepressure = porepressure
  coupling_type = Hydro
  gravity = '0 0 0'
  fp = the_simple_fluid
  mass_fraction_vars = tracer
[]

[BCs]
  [./constant_injection_porepressure]
    type = PresetBC
    variable = porepressure
    value = 1
    boundary = left
  [../]
  [./no_tracer_on_left]
    type = PresetBC
    variable = tracer
    value = 0
    boundary = left
  [../]
  [./remove_component_1]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure
    boundary = right
    fluid_phase = 0
    pt_vals = '0 1E3'
    multipliers = '0 1E3'
    mass_fraction_component = 1
    use_mobility = true
    flux_function = 1E3
  [../]
  [./remove_component_0]
    type = PorousFlowPiecewiseLinearSink
    variable = tracer
    boundary = right
    fluid_phase = 0
    pt_vals = '0 1E3'
    multipliers = '0 1E3'
    mass_fraction_component = 0
    use_mobility = true
    flux_function = 1E3
  [../]
[]

[Modules]
  [./FluidProperties]
    [./the_simple_fluid]
      type = SimpleFluidProperties
      bulk_modulus = 2E9
      viscosity = 1.0
      density0 = 1000.0
    [../]
  [../]
[]

[Materials]
  [./porosity_qp]
    type = PorousFlowPorosity
    porosity_zero = 0.1
  [../]
  [./porosity]
    type = PorousFlowPorosity
    porosity_zero = 0.1
    at_nodes = true
  [../]
  [./permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1E-2 0 0   0 1E-2 0   0 0 1E-2'
  [../]
[]

[Preconditioning]
  active = basic
  [./basic]
    type = SMP
    full = true
    petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = ' asm      lu           NONZERO                   2'
  [../]
  [./preferred_but_might_not_be_installed]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  [../]
[]

[VectorPostprocessors]
  [./tracer]
    type = LineValueSampler
    start_point = '0 0 0'
    end_point = '1 0 0'
    num_points = 101
    sort_by = x
    variable = tracer
  [../]
[]

[Executioner]
  type = Transient
  solve_type = Newton
  end_time = 6
  dt = 6E-1
  nl_abs_tol = 1E-8
  timestep_tolerance = 1E-3
[]

[Outputs]
  #exodus = true
  csv = true
  execute_on = final
[]

