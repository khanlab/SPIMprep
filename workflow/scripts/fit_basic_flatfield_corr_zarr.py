from random import shuffle
from basicpy import BaSiC
from skimage.transform import resize
import dask.array as da
import zarr
from dask.diagnostics import ProgressBar

in_zarr = snakemake.input.zarr
channel=int(snakemake.params.channel)
max_n_images=int(snakemake.params.max_n_images)

arr = da.from_zarr(in_zarr)
shape = arr.shape
n_tiles=shape[0]
n_chans=shape[1]
n_slices=shape[2]
n_x=shape[3]
n_y=shape[4]

#select the channel first
arr = arr[:,channel,:,:,:]

arr = da.reshape(arr,(n_tiles*n_slices,n_x,n_y))


#if number of images less than the max, then use entire array
if arr.shape[0] > max_n_images:
    #want to pick a random sample of max_n_images indices
    indices = [i for i in range(arr.shape[0])]
    shuffle(indices)
    arr = arr[indices[:max_n_images],:,:]

#now we want to apply resizing to all the images in the array
#but do it in parallel

def resize_parallel(x):
    return resize(x.squeeze(),(128,128),preserve_range=True).reshape((1,128,128))

#since each slice is a chunk/block, can use map_blocks to apply the function to all blocks
arr_resized = da.map_blocks(resize_parallel,arr)
   
print('concatenating and resizing with dask')
with ProgressBar():
    concat_images = arr_resized.compute()

basic = BaSiC(max_workers=snakemake.threads,
                **snakemake.params.basic_opts
                )


print('fitting BaSiC...')
basic.fit(concat_images)
print('fitting BaSiC...done')


#add plot and save model here
print('saving model')
basic.save_model(snakemake.output.model_dir)




