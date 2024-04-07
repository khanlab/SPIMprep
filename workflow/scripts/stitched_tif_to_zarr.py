from glob import glob
import dask.array as da
import dask.array.image 
from itertools import product
from dask.diagnostics import ProgressBar

rechunk_size=snakemake.params.rechunk_size


in_tif_glob = snakemake.params.in_tif_glob.format(stain=snakemake.wildcards.stain)


darr = dask.array.image.imread(in_tif_glob)
        

#rescale intensities, and recast  -- only recast
#darr = darr * snakemake.params.intensity_rescaling
darr = darr.astype('uint16')

#now we can do the computation itself, storing to zarr
print('writing images to zarr with dask')
with ProgressBar():
    da.to_zarr(darr,snakemake.output.zarr,component='fused/s0',overwrite=True,dimension_separator='/')


