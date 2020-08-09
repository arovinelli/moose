[Mesh]
  [./gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 3
    ny = 4
    nz = 5
    xmin = -1
    xmax = 1
    ymin = -2
    ymax = 2
    zmin = -3
    zmax = 3
  []

  [./subdomain0]
    type = AddSideSetsForPeriodicBCGenerator
    input = gmg
    dx_dy_dz = '2 4 6'
  []
[]

[Outputs]
  exodus = true
[]
