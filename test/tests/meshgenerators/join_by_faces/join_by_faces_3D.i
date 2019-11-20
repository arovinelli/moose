[Mesh]
  [./msh]
    type = FileMeshGenerator
    file = break_mesh_3D_auto_in.e
  []
  [./join]
    input = msh
    type = JoinMeshByFaces
  []
[]

[Outputs]
  exodus = true
[]
