[Mesh]
  ## This file is where you split your mesh
  ## execution must NOT be --mes-only as we are saving element integers for later use
  [msh]
    type = FileMeshGenerator
    file = RVE_test_mesh.e
  []
  [add_side_sets]
    input = msh
    type = SideSetsFromNormalsGenerator
    normals = '0 -1  0
               0  1  0
               -1 0  0
               1  0  0
               0  0 -1
               0  0  1'
    fixed_normal = true
    new_boundary = 'y0 y1 x0 x1 z0 z1'
  []
[]

[Problem]
  kernel_coverage_check = false
  solve = false
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
[]
