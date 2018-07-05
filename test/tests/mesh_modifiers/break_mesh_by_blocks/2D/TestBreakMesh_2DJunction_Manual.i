[Mesh]
  file = 4ElementJunction.e
[]

[MeshModifiers]
  [./breakmesh]
    type = BreakMeshByBlockManual_2DJunction
  [../]
[]

# This input file is intended to be run with the "--mesh-only" option so
# no other sections are required
