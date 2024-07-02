from jinja2 import Environment, FileSystemLoader
import os
from ome_zarr.io import parse_url
from ome_zarr.reader import Reader
import matplotlib.pyplot as plt
import numpy as np
import math

# load the html template from jinja
file_loader = FileSystemLoader(".")
env = Environment(loader=file_loader)
template = env.get_template("qc/resources/ff_html_temp.html")

# User set configurations
ff_s_start=snakemake.params.ff_s_start
ff_s_step=snakemake.params.ff_s_step
ff_cmap=snakemake.params.ff_cmap

# input zarr files
ff_corr = snakemake.input.corr
ff_uncorr = snakemake.input.uncorr

# output files
out_html = snakemake.output.html
corr_image_dir = snakemake.output.corr_images_dir
uncorr_image_dir = snakemake.output.uncorr_images_dir

# Read the corrected and uncorrected ome-zarr data
proc_reader= Reader(parse_url(ff_corr))
unproc_reader = Reader(parse_url(ff_uncorr))
proc_data=list(proc_reader())[0].data
unproc_data = list(unproc_reader())[0].data

# create directories for corrected and uncorrected images
os.makedirs(corr_image_dir, exist_ok=True)
os.mkdir(uncorr_image_dir)


chunks = []

# Create a list to store slices ordered by channel and tile
for chunk,(tile_corr, tile_uncorr) in enumerate(zip(proc_data[0], unproc_data[0])):
    channels = []
    for chan_num, (channel_corr, channel_uncorr) in enumerate(zip(tile_corr, tile_uncorr)):
        slice = ff_s_start
        slices = []
        while(slice<len(channel_corr)):
            # Get the contrast limits by removing the absolute highest and lowest data
            sorted_array = np.sort(np.array(channel_corr[slice].flatten()))[::-1]
            cmax = sorted_array[math.floor(len(sorted_array)*1/100)]
            cmin = sorted_array[math.floor(len(sorted_array)*99/100)]
            
            # clip and plot data then save image to wanted file path
            clipped_data_corr = np.clip(channel_corr[slice],cmin,cmax)
            clipped_data_uncorr = np.clip(channel_uncorr[slice], cmin, cmax)

            # get image names and paths
            corr_image_name = f"chunk-{chunk}-channel-{chan_num}-slice-{slice}.jpg"
            corrected_img_path = f"images/corr/{corr_image_name}"
            uncorr_image_name = f"chunk-{chunk}-channel-{chan_num}-slice-{slice}.jpg"
            uncorrected_img_path = f"images/uncorr/{uncorr_image_name}"

            # Save images 
            plt.imsave(corr_image_dir+"/"+corr_image_name, clipped_data_corr, cmap=ff_cmap)
            plt.imsave(uncorr_image_dir+"/"+uncorr_image_name, clipped_data_uncorr, cmap=ff_cmap)

            # create an object to store key image info
            image = {"slice": slice, "img_corr": corrected_img_path, "img_uncorr": uncorrected_img_path}
            slices.append(image)
            slice+=ff_s_step
        channels.append(slices)
    chunks.append(channels)

# pass the chunks array to the template to render the html
output = template.render(chunks=chunks, numColumns=3)
# Write out html file
with open(out_html, 'w') as f:
    f.write(output)
            


