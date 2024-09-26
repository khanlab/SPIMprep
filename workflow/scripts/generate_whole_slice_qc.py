import os
import math
import dask.array as da
import matplotlib.pyplot as plt
import numpy as np
from jinja2 import Environment, FileSystemLoader
from upath import UPath as Path
from lib.cloud_io import get_fsspec, is_remote
import zarr

# load jinja html template
file_loader = FileSystemLoader(".")
env = Environment(loader=file_loader)
template = env.get_template(snakemake.input.ws_html)

# user set configurations
ws_s_step=snakemake.params.ws_s_step
ws_s_start=snakemake.params.ws_s_start
ws_cmap=snakemake.params.ws_cmap

# output paths
image_dir = snakemake.output.images_dir
out_html = snakemake.output.html


uri = snakemake.params.uri
if is_remote(uri):
    fs_args={'storage_provider_settings':snakemake.params.storage_provider_settings,'creds':snakemake.input.creds}
else:
    fs_args={}

fs = get_fsspec(uri,**fs_args)

if Path(uri).suffix == '.zip':
    store = zarr.storage.ZipStore(Path(uri).path,dimension_separator='/',mode='r')
else:
    store = zarr.storage.FSStore(Path(uri).path,fs=fs,dimension_separator='/',mode='r')

darr = da.from_zarr(store,component='/5')

os.makedirs(image_dir, exist_ok=True)

channels = []
# for each channel add the images of the most downsampled data
for chan_num, channel in enumerate(darr):
    slice = ws_s_start
    chan = []
    slices = []
    while(slice<len(channel)):
        # Get contrast limits
        sorted_array = np.sort(np.array(channel[slice].flatten()))[::-1]
        cmax = sorted_array[math.floor(len(sorted_array)*1/100)]
        cmin = sorted_array[math.floor(len(sorted_array)*99/100)]

        # normalize data and save image 
        normalize_data= np.clip(channel[slice],cmin,cmax)
        img_name = f"channel-{chan_num}-slice-{slice}.jpg"
        proc_img_path = f"images/whole/{img_name}"
        plt.imsave(image_dir+"/"+img_name, normalize_data, cmap=ws_cmap)

        # create object containing slice image info
        slices.append({"slice": slice, "img_path": proc_img_path})
        slice += ws_s_step
        # every list here represents a line of images
        if(len(slices) == 5):
            chan.append(slices)
            slices = []
        elif slice >=len(channel):
            chan.append(slices)
    # full list containing all slices ordered by channel
    channels.append(chan)

# render the template and write out file
output = template.render(channels = channels)
with open(out_html, "w") as f:
    f.write(output)
