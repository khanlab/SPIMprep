import h5py
import hdf5plugin
import zarr
from sys import stdout #change this to log file later..

source = h5py.File(snakemake.input.ims, mode='r')

store = zarr.ZipStore(snakemake.output.zarr,dimension_separator='/',mode='x') 
dest = zarr.group(store)

zarr.copy(source['DataSet/ResolutionLevel 0/TimePoint 0/Channel {chan}/Data'.format(chan=snakemake.params.channel)], dest,log=stdout,compressor=None)
source.close()

