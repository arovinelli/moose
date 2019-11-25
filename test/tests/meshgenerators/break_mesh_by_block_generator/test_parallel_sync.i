[Mesh]
  [./fmg]
    type = FileMeshGenerator
    file = coh3D_3Blocks.e
  []

  [./breakmesh]
    type = ParallelSyncExample
    input = fmg
  []
[]

[Outputs]
  exodus = true
[]
