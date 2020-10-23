#!/bin/bash

NCPU="8"
CHECKPOINT_FILE_NAME="testfile.cpa"
MOOSE_TEST_EXEXCUTIONER_TYPE="dbg"


# ## this is tthe first step taht normally isn't needed, it simply generates a single exodus fiel that will be used by step 2
# $MOOSE_DIR/test/moose_test-$MOOSE_TEST_EXEXCUTIONER_TYPE -i simple_diffusion_test_mesh.i --mesh-only
#
## Step 2: here we are breaking the mesh in SERIAL and saving info we will need later
$MOOSE_DIR/test/moose_test-$MOOSE_TEST_EXEXCUTIONER_TYPE -i simple_diffusion_test_mesh_break.i --split-mesh $NCPU --split-file $CHECKPOINT_FILE_NAME

# ## Step 3: actual mesh split, the generated checkpoint files  will be used to run the simualtion
# $MOOSE_DIR/test/moose_test-$MOOSE_TEST_EXEXCUTIONER_TYPE -i simple_diffusion_test_mesh_prepare_split.i Mesh/msh/file=simple_diffusion_test_mesh_break_out.e Mesh/msh/exodus_extra_element_integers='bmbb_element_id' --split-mesh $NCPU --split-file $CHECKPOINT_FILE_NAME
# #

##Step 4: runnignt the simulation
mpirun -np $NCPU $MOOSE_DIR/test/moose_test-$MOOSE_TEST_EXEXCUTIONER_TYPE -i simple_diffusion_test_mesh_prepare_split.i --distributed-mesh Mesh/msh/file=$CHECKPOINT_FILE_NAME Outputs/nemesis=true Outputs/checkpoint=true

##Step 5: restarting the simulation using the msae input file, notice the problem paramter to relaod stateful var
RESTART_CHECKPOINT_FILE_LOCATION="simple_diffusion_test_mesh_prepare_split_out_cp"
CHECKPOINT_ID="0001"
mpirun -np $NCPU $MOOSE_DIR/test/moose_test-$MOOSE_TEST_EXEXCUTIONER_TYPE -i simple_diffusion_test_mesh_prepare_split.i --distributed-mesh Mesh/msh/file=$RESTART_CHECKPOINT_FILE_LOCATION/${CHECKPOINT_ID}_mesh.cpr  Problem/restart_file_base=$RESTART_CHECKPOINT_FILE_LOCATION/$CHECKPOINT_ID Outputs/nemesis=true Outputs/checkpoint=true Outputs/file_base="restart"

##Step 5b: restartwith same output file name: DON'T DO THIS, data in teh exodus file will be lost
RESTART_CHECKPOINT_FILE_LOCATION="simple_diffusion_test_mesh_prepare_split_out_cp"
CHECKPOINT_ID="0001"
mpirun -np $NCPU $MOOSE_DIR/test/moose_test-$MOOSE_TEST_EXEXCUTIONER_TYPE -i simple_diffusion_test_mesh_prepare_split.i --distributed-mesh Mesh/msh/file=$RESTART_CHECKPOINT_FILE_LOCATION/${CHECKPOINT_ID}_mesh.cpr  Problem/restart_file_base=$RESTART_CHECKPOINT_FILE_LOCATION/$CHECKPOINT_ID Outputs/nemesis=true Outputs/checkpoint=true

##Step 5, comaprison with monolitic mesh
mpirun -np $NCPU $MOOSE_DIR/test/moose_test-$MOOSE_TEST_EXEXCUTIONER_TYPE -i simple_diffusion_monolithic.i
