import tifffile
import xmltodict 
import json
import re
from itertools import product
from snakemake.io import glob_wildcards

in_tif_pattern = snakemake.params.in_tif_pattern

#add a wildcard constraint to ensure no
#subfolders get parsed (ie don't match anything with / in it):
prefix_constraint=r'[^/]+' 
in_tif_pattern_constrained = in_tif_pattern.replace('{prefix}',f'{{prefix,{prefix_constraint}}}')

#parse the filenames to get number of channels, tiles etc..
prefix, tilex, tiley, channel, zslice = glob_wildcards(in_tif_pattern_constrained)

tiles_x = sorted(list(set(tilex)))
tiles_y = sorted(list(set(tiley)))
channels = sorted(list(set(channel)))
zslices = sorted(list(set(zslice)))
prefixes = sorted(list(set(prefix)))

#read in series metadata from first file
in_tif = in_tif_pattern.format(tilex=tiles_x[0],tiley=tiles_y[0],prefix=prefixes[0],channel=channels[0],zslice=zslices[0])

raw_tif = tifffile.TiffFile(in_tif,mode='r')

axes = raw_tif.series[0].get_axes()
shape = raw_tif.series[0].get_shape()

#axes normally should be CZYX - if not, then we fail out
if not axes == 'CZYX':
    print('ERROR: cannot deal with axes other than CZYX')
    #sys.exit(1)

ome_dict = xmltodict.parse(raw_tif.ome_metadata)
physical_size_x = ome_dict['OME']['Image']['Pixels']['@PhysicalSizeX']
physical_size_y = ome_dict['OME']['Image']['Pixels']['@PhysicalSizeY']
physical_size_z = ome_dict['OME']['Image']['Pixels']['@PhysicalSizeZ']
custom_metadata = ome_dict['OME']['Image']['ca:CustomAttributes']



#read tile configuration from the microscope metadata

tile_config_pattern=r"Blaze\[(?P<tilex>[0-9]+) x (?P<tiley>[0-9]+)\]_C(?P<channel>[0-9]+)_xyz-Table Z(?P<zslice>[0-9]+).ome.tif;;\((?P<x>\S+), (?P<y>\S+),(?P<chan>\S+), (?P<z>\S+)\)"
tile_pattern = re.compile(tile_config_pattern)

#put it in 3 maps, one for each coord, indexed by tilex, tiley, channel, and aslice
map_x=dict()
map_y=dict()
map_z=dict()

map_tiles_to_chunk=dict()
chunks = []
for chunk,(tilex,tiley) in enumerate(product(tiles_x,tiles_y)):
    map_tiles_to_chunk[tilex+tiley] = chunk
    chunks.append(chunk)

for line in  custom_metadata['TileConfiguration']['@TileConfiguration'].split('  ')[1:]:
    
    d = re.search(tile_pattern,line).groupdict()
    chunk = map_tiles_to_chunk[d['tilex']+d['tiley']] # want the key to have chunk instad of tilex,tiley, so map to that first
    
    #key is:  tile-{chunk}_chan-{channel}_z-{zslice} 
    key = f"tile-{chunk}_chan-{d['channel']}_z-{d['zslice']}" 

    map_x[key] = float(d['x'])
    map_y[key] = float(d['y'])
    map_z[key] = float(d['z'])
    

metadata={}
metadata['tiles_x'] = tiles_x
metadata['tiles_y'] = tiles_y
metadata['channels'] = channels
metadata['zslices'] = zslices
metadata['prefixes'] = prefixes
metadata['chunks'] = chunks
metadata['axes'] = axes
metadata['shape'] = shape
metadata['physical_size_x'] = float(physical_size_x)
metadata['physical_size_y'] = float(physical_size_y)
metadata['physical_size_z'] = float(physical_size_z)
metadata['lookup_tile_offset_x'] = map_x
metadata['lookup_tile_offset_y'] = map_y
metadata['lookup_tile_offset_z'] = map_z
metadata['ome_full_metadata'] = ome_dict
metadata['PixelSize'] = [ metadata['physical_size_z']/1000.0, metadata['physical_size_y']/1000.0, metadata['physical_size_x']/1000.0 ] #zyx since OME-Zarr is ZYX
metadata['PixelSizeUnits'] = 'mm' 

#write metadata to json
with open(snakemake.output.metadata_json, 'w') as fp:
    json.dump(metadata, fp,indent=4)
