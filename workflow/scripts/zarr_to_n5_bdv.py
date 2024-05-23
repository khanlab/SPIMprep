import json
import numpy as np
import npy2bdv
import tifffile
import zarr
import dask
import dask.array as da
from pathlib import Path
from shutil import rmtree
from itertools import product
from dask.diagnostics import ProgressBar
import xml.etree.ElementTree as ET

def update_xml_h5_to_n5(in_xml,out_xml,in_n5):
    """updates the dataset.xml file 
    in-place to replace hdf5 with n5"""
    tree = ET.parse(in_xml)

    root = tree.getroot()
    for elem in root.iter():
        if 'format' in elem.keys():
            elem.set('format','bdv.n5')
            elem.set('version','1.0')
            for child_elem in elem.iter('hdf5'):
                child_elem.tag = 'n5'
                child_elem.text = str(Path(in_n5).name)

    tree.write(out_xml,encoding='utf-8')
    return


in_zarr = snakemake.input.zarr

max_downsampling_layers=snakemake.params.max_downsampling_layers

#load data  (tiles,chans,zslices,x,y)
darr = da.from_zarr(in_zarr)

(n_tiles,n_chans,n_z,n_x,n_y) = darr.shape

#read metadata json
with open(snakemake.input.metadata_json) as fp:
    metadata = json.load(fp)


#create temp_folder for h5/xml
temp_bdv_dir=Path(snakemake.params.temp_h5).parent
temp_bdv_dir.mkdir(parents=True, exist_ok=True)

#create new bdv volume
bdv_writer = npy2bdv.BdvWriter(snakemake.params.temp_h5,
                                overwrite=True,
                                nchannels=len(metadata['channels']), 
                                ntiles=len(metadata['tiles_x'])*len(metadata['tiles_y']),
                                blockdim=((1,256,256),))

bdv_writer.set_attribute_labels('channel', metadata['channels'])

    
for i_chan,channel in enumerate(metadata['channels']):
    for i_tile in metadata['chunks']:
        
        key = f'tile-{i_tile}_chan-{channel}_z-0000'
        affine = np.eye(3,4)
        affine[0,3] = metadata['lookup_tile_offset_x'][key] / float(metadata['physical_size_x'])
        affine[1,3] = -metadata['lookup_tile_offset_y'][key] / float(metadata['physical_size_y'])
        affine[2,3] = -metadata['lookup_tile_offset_z'][key] / float(metadata['physical_size_z'])
        bdv_writer.append_view(stack=None, 
                                virtual_stack_dim=(n_z,n_y,n_x),
                                tile=i_tile,channel=i_chan,
                                voxel_units='μm',
                                voxel_size_xyz=(metadata['physical_size_x'],
                                               metadata['physical_size_y'],
                                               metadata['physical_size_z']),
                                m_affine=affine,
                                name_affine='Translation to Regular Grid',
                                calibration=(1,1,float(metadata['physical_size_z'])/float(metadata['physical_size_x'])),
                              )




print('writing dataset xml and (empty h5) using npy2bdv')
bdv_writer.write_xml()
bdv_writer.close()

print('updating xml to point to n5 instead of h5')
update_xml_h5_to_n5(snakemake.params.temp_xml,snakemake.output.bdv_xml,snakemake.output.bdv_n5)

print('removing empty bdv h5/xml')
rmtree(temp_bdv_dir)

print('writing data to n5')
n5_store = zarr.n5.N5Store(snakemake.output.bdv_n5)

for setup_i,((chan_i,chan),(tile_i,tile)) in enumerate(product(enumerate(metadata["channels"]),enumerate(metadata["chunks"]))):

    ds_list=[] #for setup-level attrs
    for ds in range(max_downsampling_layers):
        step=2**ds  #1,2,4,8.. 
        zstack = da.squeeze(darr[tile_i,chan_i,:,::step,::step])
        print(f'writing to setup{setup_i}/timepoint0/s{ds}')
        with ProgressBar():
            zstack.to_zarr(n5_store,component=f'setup{setup_i}/timepoint0/s{ds}',overwrite=True,compute=True) 
        #add attributes for downsampling
        a = zarr.open_array(store=n5_store,path=f'setup{setup_i}/timepoint0/s{ds}',mode='r+')
        a.attrs['downsamplingFactors']=[step,step,1]
        ds_list.append([step,step,1])

    #add attributes for downsampling as a list, and datatype to the setup# level
    g = zarr.open_group(store=n5_store,path=f'setup{setup_i}',mode='r+')
    g.attrs['downsamplingFactors']=ds_list
    g.attrs['dataType']='uint16'


