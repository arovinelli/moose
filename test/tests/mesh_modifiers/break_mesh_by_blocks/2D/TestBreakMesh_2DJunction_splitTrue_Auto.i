[Mesh]
  file = 4ElementJunction.e
[]

[MeshModifiers]
  [./breakmesh]
    type = BreakMeshByBlock
    split_interface = true
  [../]
[]

# This input file is intended to be run with the "--mesh-only" option so
# no other sections are required
