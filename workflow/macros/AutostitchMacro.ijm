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
