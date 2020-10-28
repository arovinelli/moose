# Exception testing of PorousFlowHysteresisOrder
# Incorrect:     liquid_phase = 1
[Mesh]
  type = GeneratedMesh
  dim = 1
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [./pp]
  [../]
[]

[PorousFlowBasicTHM]
  porepressure = pp
  fp = simple_fluid
[]

[Modules]
  [./FluidProperties]
    [./simple_fluid]
      type = SimpleFluidProperties
    [../]
  [../]
[]


[Materials]
  [./porosity]
    type = PorousFlowPorosity
    porosity_zero = 0.1
  [../]
  [./biot_modulus]
    type = PorousFlowConstantBiotModulus
    biot_coefficient = 0.8
    solid_bulk_compliance = 2e-7
    fluid_bulk_modulus = 1e7
  [../]
  [./permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1e-13 0 0   0 1e-13 0   0 0 1e-13'
  [../]
  [./hys_order]
    type = PorousFlowHysteresisOrder
    liquid_phase = 1
  [../]
[]

[Preconditioning]
  [./basic]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = Newton
  dt = 1
  end_time = 1
[]
