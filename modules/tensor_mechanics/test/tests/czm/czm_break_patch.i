# Patch test for cohesive zone modeling to check convergence
[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [./msh]
    type = FileMeshGenerator
    file = patch_mesh.e
  []
  [./transform]
    type = TransformGenerator
    input = msh
    transform = TRANSLATE
    vector_value = '-0.5 -0.5 -0.5'
  []
  [./split]
    type = BreakMeshByBlockGenerator
    input = transform
    write_fake_neighbor_list_to_file = true
  []
[]

[Outputs]
  exodus = true
[]
