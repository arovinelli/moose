CDF      
      
len_string     !   len_line   Q   four      	time_step          len_name   !   num_dim       	num_nodes      	   num_elem      
num_el_blk        num_node_sets         num_side_sets         num_el_in_blk1        num_nod_per_el1       num_el_in_blk2        num_nod_per_el2       num_side_ss1      num_side_ss2      num_side_ss3      num_side_ss4      num_side_ss5      num_nod_ns1       num_nod_ns2       num_nod_ns3       num_nod_ns4       num_nod_ns5       num_nod_var       num_glo_var       num_info  j         api_version       @�
=   version       @�
=   floating_point_word_size            	file_size               int64_status             title         interface_user_object_out.e    maximum_name_length                 %   
time_whole                            �(   	eb_status                             
l   eb_prop1               name      ID              
t   	ns_status         	                    
|   ns_prop1      	         name      ID              
�   	ss_status         
                    
�   ss_prop1      
         name      ID              
�   coordx                      H      
�   coordy                      H         eb_names                       D      \   ns_names      	                 �      �   ss_names      
                 �      H   
coor_names                         D      �   node_num_map                    $      4   connect1                  	elem_type         QUAD4               X   connect2                  	elem_type         QUAD4         0      h   elem_num_map                          �   elem_ss1                          �   side_ss1                          �   elem_ss2                          �   side_ss2                          �   elem_ss3                          �   side_ss3                          �   elem_ss4                          �   side_ss4                          �   elem_ss5                          �   side_ss5                          �   node_ns1                          �   node_ns2                             node_ns3                             node_ns4                             node_ns5                          (   vals_nod_var1                          H      �0   name_nod_var                       $      4   name_glo_var                       D      X   vals_glo_var                             �x   info_records                      r�      �                                                                           ��                      ��      ?�      ?�              ��      ?�      ��      ��                      ��              ?�      ?�      ?�                                                                          bottom                           left                                                              right                            top                                 bottom                           right                            top                              left                             interface                                                                                                                          	                                             	                                                                                                            	         	   u                                   mean_mat_prop                    mean_var_jump                      ####################                                                             # Created by MOOSE #                                                             ####################                                                             ### Command Line Arguments ###                                                    -i interface_user_object.i### Version Info ###                                  Framework Information:                                                           MOOSE Version:           git commit 87009251b7 on 2018-06-12                     LibMesh Version:         ceb1e02bee97b3758bfcdc37cb87bc8954a89030                PETSc Version:           3.7.6                                                   Current Time:            Tue Jun 12 13:35:10 2018                                Executable Timestamp:    Tue Jun 12 13:34:13 2018                                                                                                                                                                                                  ### Input File ###                                                                                                                                                []                                                                                 inactive                       =                                                 initial_from_file_timestep     = LATEST                                          initial_from_file_var          = INVALID                                         element_order                  = AUTO                                            order                          = AUTO                                            side_order                     = AUTO                                            type                           = GAUSS                                         []                                                                                                                                                                [BCs]                                                                                                                                                               [./all]                                                                            boundary                     = '0 1 2 3'                                         control_tags                 = INVALID                                           enable                       = 1                                                 extra_matrix_tags            = INVALID                                           extra_vector_tags            = INVALID                                           implicit                     = 1                                                 inactive                     =                                                   isObjectAction               = 1                                                 matrix_tags                  = system                                            type                         = FunctionDirichletBC                               use_displaced_mesh           = 0                                                 variable                     = u                                                 vector_tags                  = nontime                                           diag_save_in                 = INVALID                                           function                     = fn_exact                                          save_in                      = INVALID                                           seed                         = 0                                               [../]                                                                          []                                                                                                                                                                [Executioner]                                                                      inactive                       =                                                 isObjectAction                 = 1                                               type                           = Steady                                          compute_initial_residual_before_preset_bcs = 0                                   contact_line_search_allowed_lambda_cuts = 2                                      contact_line_search_ltol       = INVALID                                         control_tags                   =                                                 enable                         = 1                                               l_abs_step_tol                 = -1                                              l_max_its                      = 10000                                           l_tol                          = 1e-05                                           line_search                    = default                                         line_search_package            = petsc                                           mffd_type                      = wp                                              nl_abs_step_tol                = 1e-50                                           nl_abs_tol                     = 1e-50                                           nl_max_funcs                   = 10000                                           nl_max_its                     = 50                                              nl_rel_step_tol                = 1e-50                                           nl_rel_tol                     = 1e-08                                           no_fe_reinit                   = 0                                               petsc_options                  = INVALID                                         petsc_options_iname            = INVALID                                         petsc_options_value            = INVALID                                         restart_file_base              =                                                 solve_type                     = NEWTON                                          splitting                      = INVALID                                       []                                                                                                                                                                [Executioner]                                                                    []                                                                                                                                                                [Functions]                                                                                                                                                         [./ffn]                                                                            inactive                     =                                                   isObjectAction               = 1                                                 type                         = ParsedFunction                                    control_tags                 = Functions                                         enable                       = 1                                                 vals                         = INVALID                                           value                        = -4                                                vars                         = INVALID                                         [../]                                                                                                                                                             [./fn_exact]                                                                       inactive                     =                                                   isObjectAction               = 1                                                 type                         = ParsedFunction                                    control_tags                 = Functions                                         enable                       = 1                                                 vals                         = INVALID                                           value                        = x*x+y*y                                           vars                         = INVALID                                         [../]                                                                          []                                                                                                                                                                [Kernels]                                                                                                                                                           [./diff]                                                                           inactive                     =                                                   isObjectAction               = 1                                                 type                         = Diffusion                                         block                        = INVALID                                           control_tags                 = Kernels                                           diag_save_in                 = INVALID                                           enable                       = 1                                                 extra_matrix_tags            = INVALID                                           extra_vector_tags            = INVALID                                           implicit                     = 1                                                 matrix_tags                  = system                                            save_in                      = INVALID                                           seed                         = 0                                                 use_displaced_mesh           = 0                                                 variable                     = u                                                 vector_tags                  = nontime                                         [../]                                                                                                                                                             [./ffn]                                                                            inactive                     =                                                   isObjectAction               = 1                                                 type                         = BodyForce                                         block                        = INVALID                                           control_tags                 = Kernels                                           diag_save_in                 = INVALID                                           enable                       = 1                                                 extra_matrix_tags            = INVALID                                           extra_vector_tags            = INVALID                                           function                     = ffn                                               implicit                     = 1                                                 matrix_tags                  = system                                            postprocessor                = 1                                                 save_in                      = INVALID                                           seed                         = 0                                                 use_displaced_mesh           = 0                                                 value                        = 1                                                 variable                     = u                                                 vector_tags                  = nontime                                         [../]                                                                          []                                                                                                                                                                [Materials]                                                                                                                                                         [./stateful1]                                                                      inactive                     =                                                   isObjectAction               = 1                                                 type                         = StatefulMaterial                                  block                        = 0                                                 boundary                     = INVALID                                           compute                      = 1                                                 constant_on                  = NONE                                              control_tags                 = Materials                                         enable                       = 1                                                 implicit                     = 1                                                 initial_diffusivity          = 5                                                 output_properties            = INVALID                                           outputs                      = none                                              seed                         = 0                                                 use_displaced_mesh           = 0                                               [../]                                                                                                                                                             [./stateful2]                                                                      inactive                     =                                                   isObjectAction               = 1                                                 type                         = StatefulMaterial                                  block                        = 1                                                 boundary                     = INVALID                                           compute                      = 1                                                 constant_on                  = NONE                                              control_tags                 = Materials                                         enable                       = 1                                                 implicit                     = 1                                                 initial_diffusivity          = 2                                                 output_properties            = INVALID                                           outputs                      = none                                              seed                         = 0                                                 use_displaced_mesh           = 0                                               [../]                                                                          []                                                                                                                                                                [Mesh]                                                                             inactive                       =                                                 displacements                  = INVALID                                         block_id                       = INVALID                                         block_name                     = INVALID                                         boundary_id                    = INVALID                                         boundary_name                  = INVALID                                         construct_side_list_from_node_list = 0                                           ghosted_boundaries             = INVALID                                         ghosted_boundaries_inflation   = INVALID                                         isObjectAction                 = 1                                               second_order                   = 0                                               skip_partitioning              = 0                                               type                           = GeneratedMesh                                   uniform_refine                 = 0                                               allow_renumbering              = 1                                               bias_x                         = 1                                               bias_y                         = 1                                               bias_z                         = 1                                               centroid_partitioner_direction = INVALID                                         construct_node_list_from_side_list = 1                                           control_tags                   =                                                 dim                            = 2                                               elem_type                      = QUAD4                                           enable                         = 1                                               gauss_lobatto_grid             = 0                                               ghosting_patch_size            = INVALID                                         max_leaf_size                  = 10                                              nemesis                        = 0                                               nx                             = 2                                               ny                             = 2                                               nz                             = 1                                               parallel_type                  = DEFAULT                                         partitioner                    = default                                         patch_size                     = 40                                              patch_update_strategy          = never                                           xmax                           = 1                                               xmin                           = -1                                              ymax                           = 1                                               ymin                           = -1                                              zmax                           = 1                                               zmin                           = 0                                             []                                                                                                                                                                [Mesh]                                                                           []                                                                                                                                                                [MeshModifiers]                                                                                                                                                     [./interface]                                                                      inactive                     =                                                   isObjectAction               = 1                                                 type                         = SideSetsBetweenSubdomains                         control_tags                 = MeshModifiers                                     depends_on                   = subdomain_id                                      enable                       = 1                                                 force_prepare                = 0                                                 master_block                 = 0                                                 new_boundary                 = interface                                         paired_block                 = 1                                               [../]                                                                                                                                                             [./subdomain_id]                                                                   inactive                     =                                                   isObjectAction               = 1                                                 type                         = AssignElementSubdomainID                          control_tags                 = MeshModifiers                                     depends_on                   = INVALID                                           element_ids                  = INVALID                                           enable                       = 1                                                 force_prepare                = 0                                                 subdomain_ids                = '0 1 1 1'                                       [../]                                                                          []                                                                                                                                                                [Outputs]                                                                          append_date                    = 0                                               append_date_format             = INVALID                                         checkpoint                     = 0                                               color                          = 1                                               console                        = 1                                               controls                       = 0                                               csv                            = 0                                               dofmap                         = 0                                               execute_on                     = 'INITIAL TIMESTEP_END'                          exodus                         = 1                                               file_base                      = INVALID                                         gmv                            = 0                                               gnuplot                        = 0                                               hide                           = INVALID                                         inactive                       =                                                 interval                       = 1                                               nemesis                        = 0                                               output_if_base_contains        = INVALID                                         print_linear_residuals         = 1                                               print_mesh_changed_info        = 0                                               print_perf_log                 = 0                                               show                           = INVALID                                         solution_history               = 0                                               sync_times                     =                                                 tecplot                        = 0                                               vtk                            = 0                                               xda                            = 0                                               xdr                            = 0                                             []                                                                                                                                                                [Postprocessors]                                                                                                                                                    [./mean_mat_prop]                                                                  inactive                     =                                                   isObjectAction               = 1                                                 type                         = InterfaceUOPPS                                    allow_duplicate_execution_on_initial = 0                                         control_tags                 = Postprocessors                                    enable                       = 1                                                 execute_on                   = TIMESTEP_END                                      force_preaux                 = 0                                                 outputs                      = INVALID                                           return_var_jump              = 0                                                 use_displaced_mesh           = 0                                                 user_object                  = interface_uo                                    [../]                                                                                                                                                             [./mean_var_jump]                                                                  inactive                     =                                                   isObjectAction               = 1                                                 type                         = InterfaceUOPPS                                    allow_duplicate_execution_on_initial = 0                                         control_tags                 = Postprocessors                                    enable                       = 1                                                 execute_on                   = TIMESTEP_END                                      force_preaux                 = 0                                                 outputs                      = INVALID                                           return_var_jump              = 1                                                 use_displaced_mesh           = 0                                                 user_object                  = interface_uo                                    [../]                                                                          []                                                                                                                                                                [UserObjects]                                                                                                                                                       [./interface_uo]                                                                   inactive                     =                                                   isObjectAction               = 1                                                 type                         = InterfaceUO                                       allow_duplicate_execution_on_initial = 0                                         boundary                     = interface                                         control_tags                 = UserObjects                                       diffusivity                  = diffusivity                                       enable                       = 1                                                 execute_on                   = 'INITIAL TIMESTEP_END'                            implicit                     = 1                                                 use_displaced_mesh           = 0                                                 use_old_prop                 = 0                                                 variable                     = u                                               [../]                                                                          []                                                                                                                                                                [Variables]                                                                                                                                                         [./u]                                                                              block                        = INVALID                                           eigen                        = 0                                                 family                       = LAGRANGE                                          inactive                     =                                                   initial_condition            = INVALID                                           order                        = FIRST                                             outputs                      = INVALID                                           scaling                      = 1                                                 initial_from_file_timestep   = LATEST                                            initial_from_file_var        = INVALID                                         [../]                                                                          []                                                                                                                                                                                 ?�      @       ?�      <�������?�      @       ?�      ?�      @       @       @              