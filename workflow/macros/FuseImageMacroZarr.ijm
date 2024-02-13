args = getArgument()
args = split(args, " ");

in_dataset_xml = args[0]
downsampling = args[1]
channel = args[2]
output_zarr = args[3]
block_size_x = args[4]
block_size_y = args[5]
block_size_z = args[6]
block_size_factor_x = args[7]
block_size_factor_y = args[8]
block_size_factor_z = args[9]

run("Fuse dataset ...",
 "select=" + in_dataset_xml +
    " process_angle=[All angles] " + 
    " process_channel=[Single channel (Select from List)] " + 
    " processing_channel=[channel " + channel + "]" +
    " process_illumination=[All illuminations] " + 
    " process_tile=[All tiles] " + 
    " process_timepoint=[All Timepoints] " + 
    " bounding_box=[Currently Selected Views] " + 
    " preserve_original " +
    " downsampling=" + downsampling + 
    " interpolation=[Linear Interpolation] " + 
    " pixel_type=[16-bit unsigned integer] " + 
    " interest_points_for_non_rigid=[-= Disable Non-Rigid =-] " + 
    " blend produce=[Each timepoint & channel] " + 
    " fused_image=[ZARR/N5/HDF5 export using N5-API] " +
    " define_input=[Auto-load from input data (values shown below)] " +
    " export=ZARR create_0 " +
    " zarr_dataset_path=" + output_zarr + 
    " zarr_base_dataset=/ " +
    " zarr_dataset_extension=/s0 " + 
    " show_advanced_block_size_options " +
    " block_size_x=" + block_size_x +
    " block_size_y=" + block_size_y +
    " block_size_z=" + block_size_z +
    " block_size_factor_x=" + block_size_factor_x +
    " block_size_factor_y=" + block_size_factor_y +
    " block_size_factor_z=" + block_size_factor_z +
    " subsampling_factors=[{ {1,1,1}}]"); 


// quit after we are finished
eval("script", "System.exit(0);");
