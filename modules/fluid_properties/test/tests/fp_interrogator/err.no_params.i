[FluidPropertiesInterrogator]
  fp = fp
[]

[Modules]
  [./FluidProperties]
    [./fp]
      type = IdealGasFluidProperties
      gamma = 1.4
      R = 286.7
      mu = 1.823e-05
      k = 0.02568
    [../]
  [../]
[]
