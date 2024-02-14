args = getArgument()
args = split(args, " ");

dataset_xml = args[0]
ds_x = args[1]
ds_y = args[2]
ds_z = args[3]
min_r = args[4]

run("Calculate pairwise shifts ...",
 "select=" + dataset_xml + 
    " process_angle=[All angles] " + 
    " process_channel=[All channels] " + 
    " process_illumination=[All illuminations] " + 
    " process_tile=[All tiles] " + 
    " process_timepoint=[All Timepoints] " + 
    " method=[Phase Correlation] " + 
    " channels=[Average Channels] " + 
    " downsample_in_x=" +  ds_x + 
    " downsample_in_y=" +  ds_y + 
    " downsample_in_z=" +  ds_z );

run("Filter pairwise shifts ...",
 "select=" + dataset_xml + 
    " min_r=" + min_r + 
    " max_r=1 " + 
    " max_shift_in_x=0 " + 
    " max_shift_in_y=0 " + 
    " max_shift_in_z=0 " + 
    " max_displacement=0");

run("Calculate pairwise shifts ...", 
"browse=" + dataset_xml + 
" select="+ dataset_xml + 
" process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] method=[Interest-Point Registration (with new Interest Points)] show_expert_grouping_options show_expert_algorithm_parameters how_to_treat_timepoints=[treat individually] how_to_treat_channels=group how_to_treat_illuminations=group how_to_treat_angles=[treat individually] how_to_treat_tiles=compare type_of_interest_point_detection=Difference-of-Gaussian label_interest_points=beads subpixel_localization=[3-dimensional quadratic fit] interest_point_specification=[Comparable to Sample & small (beads)]" + 
 " downsample_xy=" + ds_x + "x " +
 " downsample_z=" + ds_z + "x " +
 " compute_on=[CPU (Java)] registration_algorithm=[Fast descriptor-based (rotation invariant)] registration_in_between_views=[Only compare overlapping views (according to current transformations)] interest_point_inclusion=[Compare all interest point of overlapping views] interest_points=beads group_channels redundancy=0 significance=10 allowed_error_for_ransac=5 number_of_ransac_iterations=Normal");

run("ICP Refinement ...", "select="+dataset_xml+ 
 " process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] icp_refinement_type=[Simple (tile registration)] global_optimization_strategy=[Two-Round: Handle unconnected tiles, remove wrong links RELAXED (5.0x / 7.0px)] downsampling=[Downsampling 4/4/2] interest=[Average Threshold] icp_max_error=[Normal Adjustment (<5px)]");


 
run("Optimize globally and apply shifts ...",
 "select=" + dataset_xml + 
    " process_angle=[All angles] " + 
    " process_channel=[All channels] " + 
    " process_illumination=[All illuminations] " + 
    " process_tile=[All tiles] " + 
    " process_timepoint=[All Timepoints] " + 
    " relative=2.500 " + 
    " absolute=3.500 " + 
    " global_optimization_strategy=[Two-Round using Metadata to align unconnected Tiles and iterative dropping of bad links] " + 
    " fix_group_0-0,");

// quit after we are finished
eval("script", "System.exit(0);");
