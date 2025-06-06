import json
from basicpy import BaSiC
from pathlib import Path
import numpy as np
from skimage.transform import resize
import dask.array as da
from dask.diagnostics import ProgressBar



in_zarr = snakemake.input.zarr

#load data
arr = da.from_zarr(in_zarr)

img_shape = (arr.shape[3],arr.shape[4])

chan_arr_list=[]

#loops over channels (one model for each channel):
for i_channel,model in enumerate(snakemake.input.model_dirs):


    #load model (doing it here instead of using load_model function as it was buggy-- fix that later)..
    model_path = Path(model)
    with open(model_path / 'settings.json' ) as fp:
        model = json.load(fp)

    profiles = np.load(model_path / 'profiles.npz')
    model["flatfield"] = profiles['flatfield']
    model["darkfield"] = profiles['darkfield']

    basic = BaSiC()
    basic.flatfield = resize(model['flatfield'],img_shape,preserve_range=True)
    basic.darkfield = resize(model['darkfield'],img_shape,preserve_range=True)
    basic.baseline = 0

    #select the channel first
    arr_chan = arr[:,i_channel,:,:,:]

    #now we want to apply correction to all images
    #define a function to map
    def apply_basic_parallel(x):
        return np.reshape(basic.transform(x.squeeze()),(1,1,img_shape[0],img_shape[1])).astype('uint16')
    arr_corr = da.map_blocks(apply_basic_parallel,arr_chan)

    chan_arr_list.append(arr_corr)

#stack along chans
arr_stacked = da.stack(chan_arr_list,axis=1).rechunk([1,1] + snakemake.params.out_chunks)

with ProgressBar():
    da.to_zarr(arr_stacked,snakemake.output.zarr,overwrite=True)
