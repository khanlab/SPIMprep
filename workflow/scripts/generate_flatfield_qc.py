import os
import math
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from ome_zarr.io import parse_url
from ome_zarr.reader import Reader

# script for creating flatfield before/after qc snapshots and html


slice_start=snakemake.params.ff_s_start
slice_step=snakemake.params.ff_s_step
colour=snakemake.params.ff_cmap
ff_corr = snakemake.input.corr
ff_uncorr = snakemake.input.uncorr
out_html = snakemake.output.html
images_dir = snakemake.output.images_dir

#create output dir for images
Path(images_dir).mkdir(parents=True,exist_ok=True)


# Get corrected and uncorrected file paths
corr_zarr = corrected_path
uncorr_zarr = uncorrected_path

# Read the corrected and uncorrected ome-zarr data
proc_reader= Reader(parse_url(corr_zarr))
unproc_reader = Reader(parse_url(uncorr_zarr))

proc_data=list(proc_reader())[0].data
unproc_data = list(unproc_reader())[0].data

# Create html file for flatfield corrected images
with open(out_html, 'w') as fw:
    fw.write(f'''
<!DOCTYPE html>
<head>
<title>FlatField Correction Check</title>
<link rel="stylesheet" href="style.css"/>
</head>
<body>
<a href="../index.html">Back</a>
<label class="expand-options">
Expand Images on Click
<input type="checkbox" id="expand">
<input type="number" id="expand_scale" value="2">
</label>
<table id="table">
<tbody>
  <h1>Before and After Flatfield Correction</h1>                 
''')
    
    # Add images into the table for each chunk and channel
    for chunk,(tile_corr, tile_uncorr) in enumerate(zip(proc_data[0], unproc_data[0])):
        for chan_num, (channel_corr, channel_uncorr) in enumerate(zip(tile_corr, tile_uncorr)):
            slice = slice_start
            fw.write(f"""      <tr>
    <td colspan={len(channel_corr/slice_step)*2}>
      <h2>Chunk - {chunk}  Channel - {chan_num}</h2>
    </td>
  </tr>
  <tr>""")
            # Process every wanted slice within a chunk and channe;
            while(slice<len(channel_corr)):
                # Get the contrast limits by removing the absolute highest and lowest data
                sorted_array = np.sort(np.array(channel_corr[slice].flatten()))[::-1]
                cmax = sorted_array[math.floor(len(sorted_array)*1/100)]
                cmin = sorted_array[math.floor(len(sorted_array)*99/100)]
                
                # clip and plot data then save image to wanted file path
                clipped_data_corr = np.clip(channel_corr[slice],cmin,cmax)
                clipped_data_uncorr = np.clip(channel_uncorr[slice], cmin, cmax)
                corrected_img_path = Path(images_dir) / f"chunk-{chunk}-channel-{chan_num}-slice-{slice}_corr.jpg"
                uncorrected_img_path = Path(images_dir) / f"chunk-{chunk}-channel-{chan_num}-slice-{slice}_uncorr.jpg"
                plt.imsave(corrected_img_path, clipped_data_corr, cmap=colour)
                plt.imsave(uncorrected_img_path, clipped_data_uncorr, cmap=colour) 
                corr_relpath = corrected_img_path.relative_to(Path(out_html).parent)
                uncorr_relpath = uncorrected_img_path.relative_to(Path(out_html).parent)
                
                # Add the images into the html format                  
                fw.write(f'''       
    <td>
        <img src={corr_relpath}></img>
        <h3>Corrected</h3>
        <p>Slice-{slice}</p>
    </td>
    <td>
        <img src={uncorr_relpath}></img>
        <h3>Uncorrected</h3>
        <p>Slice-{slice}</p>
    </td>''')
                # Increase by user given slice step
                slice += slice_step

            fw.write("      </tr>")

    fw.write("""
  </tr>
</tbody>
</table>
<script src=image_expand.js></script>                            
             """)
            


