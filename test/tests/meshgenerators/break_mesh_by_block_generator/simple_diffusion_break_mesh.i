[Mesh]
  [msh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 2
  []
  [./subdomain_1]
    type = SubdomainBoundingBoxGenerator
    input = msh
    bottom_left = '0 0.5 0'
    top_right = '1 1 0'
    block_id = 1
  []
  [break]
    type = BreakMeshByBlockGenerator
    input = subdomain_1
    write_fake_neighbor_list_to_file = true
    fake_neighbor_list_file_name = 'simple_diffusion_presplit.txt'
  []
[]
