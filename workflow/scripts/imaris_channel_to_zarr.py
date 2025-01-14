import h5py
import hdf5plugin
import zarr
from sys import stdout #change this to log file later..

source = h5py.File(snakemake.input.ims, mode='r')
dest = zarr.open_group(snakemake.output.zarr, mode='w')
zarr.copy(source['DataSet/ResolutionLevel 0/TimePoint 0/Channel {chan}/Data'.format(chan=snakemake.params.channel)], dest,log=stdout,compressor=None)
source.close()

