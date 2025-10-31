if __name__ == "__main__":

    from dask.distributed import Client, LocalCluster

    cluster = LocalCluster(
        n_workers=int(snakemake.threads / 2),  # or 32, depending on workload
        threads_per_worker=2,  # isolate GIL
        memory_limit="auto",  # or tune to your RAM
        dashboard_address=":8788",
    )
    client = Client(cluster)
    print(cluster.dashboard_link)



    import json
    import os
    import pickle

    import dask.array as da
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    import zarr
    from multiview_stitcher import fusion, io, msi_utils, ngff_utils, registration
    from multiview_stitcher import spatial_image_utils as si_utils
    from multiview_stitcher import vis_utils
    from skimage.filters import threshold_otsu

    matplotlib.use("agg")


    def get_background_threshold(tile_means):

        from scipy.signal import find_peaks

        # Compute histogram
        hist, bin_edges = np.histogram(tile_means, bins=100)

        # Find valleys (peaks in negative histogram)
        valley_indices, _ = find_peaks(-hist)

        # Get the threshold as the bin edge corresponding to the first valley
        if valley_indices.size > 0:
            valley_threshold = bin_edges[valley_indices[0]]
        else:
            valley_threshold = None
            print("No valley detected in histogram.")

        # Plot the histogram and mark the valley threshold
        #    plt.figure(figsize=(10, 6))
        #    plt.hist(data, bins=100, edgecolor='black', alpha=0.7)
        #    if valley_threshold is not None:
        #        plt.axvline(valley_threshold, color='red', linestyle='--', label=f'Valley Threshold ≈ {valley_threshold:.2f}')
        #        plt.legend()
        #    plt.title('Histogram with Valley Threshold')
        #    plt.xlabel('Value')
        #    plt.ylabel('Frequency')
        #    plt.grid(True)
        #    plt.show()
        return valley_threshold


    in_zarr = snakemake.input.zarr

    darr = da.from_zarr(in_zarr)

    (n_tiles, n_chans, n_z, n_x, n_y) = darr.shape

    # read metadata json
    with open(snakemake.input.metadata_json) as fp:
        metadata = json.load(fp)


    msims = []
    zarr_paths = []

    curr_transform_key = "affine_metadata"
    new_transform_key = "affine_registered"


    channel = metadata["channels"][snakemake.params.reg_channel_index]

    # get tile means for calculating whether a tile contains the sample or not


    # Select the relevant channel across all tiles
    # Resulting shape: (num_tiles, z, y, x)
    selected = darr[:, snakemake.params.reg_channel_index, :, :, :]

    # Compute mean across spatial dimensions per tile
    # Resulting shape: (num_tiles,)
    print("calculating means per tile")
    tile_means = selected.mean(axis=(1, 2, 3)).compute()

    background_threshold = get_background_threshold(tile_means)
    print(f"background threshold: {background_threshold}")

    # Determine sample presence based on threshold
    contain_sample = list(tile_means > background_threshold)

    # Optional: print tile means
    for i_tile, mean_val in zip(metadata["chunks"], tile_means):
        print(f"{i_tile}, mean={mean_val}")


    for i_tile in metadata["chunks"]:

        key = f"tile-{i_tile}_chan-{channel}_z-0000"

        # read tile image
        im_data = da.squeeze(darr[i_tile, :, :, :, :])

        print(im_data.shape)
        sim = si_utils.get_sim_from_array(
            im_data,
            dims=["c", "z", "y", "x"],
            scale={
                "z": metadata["physical_size_z"],
                "y": metadata["physical_size_y"],
                "x": metadata["physical_size_x"],
            },
            translation={
                "z": -metadata["lookup_tile_offset_z"][key],
                "y": -metadata["lookup_tile_offset_y"][key],
                "x": metadata["lookup_tile_offset_x"][key],
            },
            c_coords=snakemake.params.channels,
            transform_key=io.METADATA_TRANSFORM_KEY,
        )

        # for next steps, we read things back as msim
        msim = msi_utils.get_msim_from_sim(sim)

        msims.append(msim)


    # exclude tiles containing only background from registration
    print("performing stitching registration")
    reg_result = registration.register(
        [msim for i, msim in enumerate(msims) if contain_sample[i]],
        reg_channel_index=snakemake.params.reg_channel_index,
        transform_key=curr_transform_key,
        new_transform_key=new_transform_key,
        pre_registration_pruning_method="keep_axis_aligned",
        post_registration_do_quality_filter=True,
        post_registration_quality_threshold=0.6,
        scheduler="threads",
        return_dict=True,
        **snakemake.params.registration_opts,
    )


    # Extract each affine and save to disk as a single .npz file
    affines = {}
    for imsim, msim in enumerate(msims):
        if contain_sample[imsim]:
            affine = np.array(
                msi_utils.get_transform_from_msim(msim, transform_key=new_transform_key)[0]
            )
        else:
            affine = np.eye(4)

        affines[f"tile_{imsim}"] = affine
        print(f"tile index {imsim}\n", affine)

    np.savez(snakemake.output.affines, **affines)
    print("Saved affines to registered_affines.npz")


    # Save the entire registration result dictionary
    with open(snakemake.output.reg_result_pkl, "wb") as f:
        pickle.dump(reg_result, f)

    print("Saved full registration result to Pickle")
