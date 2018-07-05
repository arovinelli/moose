# Break Mesh By Block

!syntax description /MeshModifiers/BreakMeshByBlock

This class implement a methods to split a monolithic mesh by blocks and follows the
method proposed by VP Nguyen [cite:Nguyen2014]

To split the mesh nodes belonging to a boundary between 2 blocks are duplicated, the mesh is split and a new sideset is added. The new sideset is the boundary where the Interface lives.



As an option, the interface can be split into $N$ different sidesets. $N$ is the number of adjacent block pairs. This is achieved by setting  `split_interface=true`. This is useful when modeling interfaces with different parameters.

## Example Input File Syntax



### Single interface

!listing test/tests/mesh_modifiers/break_mesh_by_blocks/2D/TestBreakMesh_2DJunction_Auto.i block=MeshModifiers

### Multiple interfaces

When `split_interface=true` the new generated interface is splitted by block pairs:

!listing test/tests/mesh_modifiers/break_mesh_by_blocks/2D/TestBreakMesh_2DJunction_splitTrue_Auto.i block=MeshModifiers

The naming convention for the interface is `bMi_bSj`, where $i$ and $j$ represents the original blocks ids.
This convention also assumes that the master interface $i$ is set on the block with lower id, while the slave interface $j$ is always on the block with higher id.
Note that the new generated interfaces might not have consecutive ids, therefore interfaces should be referred by name.

!syntax inputs /MeshModifiers/BreakMeshByBlock

!syntax parameters /MeshModifiers/BreakMeshByBlock

!bibtex bibliography
