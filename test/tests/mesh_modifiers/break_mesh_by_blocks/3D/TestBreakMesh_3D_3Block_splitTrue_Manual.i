[Mesh]
  file = coh3D_3Blocks.e
[]

[MeshModifiers]
  [./breakmesh]
    type = BreakMeshByBlockManual_3Blocks
    split_interface = true
  [../]
[]

# This input file is intended to be run with the "--mesh-only" option so
# no other sections are required
