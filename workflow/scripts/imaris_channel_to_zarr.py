from sys import stdout  # change this to log file later..

import h5py
import hdf5plugin
import zarr

source = h5py.File(snakemake.input.ims, mode="r")

store = zarr.ZipStore(snakemake.output.zarr, dimension_separator="/", mode="x")
dest = zarr.group(store)


def print_tree_item(name, obj):
    """
    A visitor function to print the name and type of HDF5 objects.
    """
    # Indentation based on group depth
    shift = name.count("/") * "  "
    if isinstance(obj, h5py.Dataset):
        print(f"{shift}Dataset: {name} (shape: {obj.shape}, dtype: {obj.dtype})")
    elif isinstance(obj, h5py.Group):
        print(f"{shift}Group: {name}/")
    else:
        print(f"{shift}Other: {name}")

    # Optionally, print attributes as well
    for key, val in obj.attrs.items():
        print(f"{shift}  Attribute: {key}: {val}")


print(
    source[
        "DataSet/ResolutionLevel 0/TimePoint 0/Channel {chan}/Data".format(
            chan=snakemake.params.channel
        )
    ].shape
)
source.visititems(print_tree_item)
print(
    source[
        "DataSet/ResolutionLevel 0/TimePoint 0/Channel {chan}/Data".format(
            chan=snakemake.params.channel
        )
    ].shape
)

# zarr.copy(
#    source[
#        "DataSet/ResolutionLevel 0/TimePoint 0/Channel {chan}/Data".format(
#            chan=snakemake.params.channel
#        )
#    ],
#    dest,
#    log=stdout,
#    compressor=None,
# )
source.close()
