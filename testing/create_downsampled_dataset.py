import numpy as np
from zarrnii import ZarrNii
import tarfile
import tifffile
import xmltodict
import dask.array as da
import os
import typer
from typing_extensions import Annotated
app = typer.Typer()

def downsample_tiff(source_dir, ds_x, ds_y, ds_z, slice_step):
    """
    Take in the original tiff data and put into zarnii to handle the downsampling
    """
    member_names = []
    members_list = []
    za = None
    with tarfile.open(source_dir, 'r') as tar:
        members = tar.getmembers()
        for member in members:
            tar.extract(member)
            member_names.append(member.name)
            members_list.append(member)
    for member in member_names:
        with tifffile.TiffFile(member) as tif:
            data = tif.series[0].asarray()
            data = da.from_array(data)
            member_slice = int(member.split("Z")[1][:4])
            channel = int(member.split("C")[1][:2]) 
            print(member_slice)
            if(za == None):
                za = ZarrNii.from_darr(data)
                meta = xmltodict.parse(tif.ome_metadata)
            elif(data.shape == (2560,2160)):
                za.darr[channel, member_slice/slice_step] = np.array(data)
            else:
                za.darr[channel, member_slice/slice_step] = np.array(data[channel][member_slice])
    za = za.downsample(along_x=ds_x, along_y=ds_y, along_z=ds_z)
    za.darr = da.from_array(np.array(za.darr).astype(np.uint16))
    return meta, za, members_list


def basic_meta_update(meta, za, ds_x=1, ds_y=1, ds_z=1):
    """
    Update the simple metadata including pixel size and the size of the array
    """
    pixel = meta['OME']['Image']['Pixels']
    pixel['@SizeX'] = f'{za.darr.shape[3]}'
    pixel['@SizeY'] = f'{za.darr.shape[2]}'
    pixel['@SizeZ'] = f'{za.darr.shape[1]}'
    pixel['@PhysicalSizeX'] = f"{float(pixel['@PhysicalSizeX'])*ds_x}"
    pixel['@PhysicalSizeY'] = f"{float(pixel['@PhysicalSizeY'])*ds_y}"
    pixel['@PhysicalSizeZ'] = f"{float(pixel['@PhysicalSizeZ'])*ds_z}"
    meta['OME']['Image']['Pixels'] = pixel
    return meta

def advanced_meta(meta, za, slice_step):
    """
    Update the tiffdata tile configuration data to ensure
    data is read and processed properly
    """
    tiff_data = meta['OME']['Image']['Pixels']['TiffData']
    new_tiff_data = []
    for single_data in tiff_data:
        slice_num = int(single_data["@FirstZ"])
        if slice_num < za.darr.shape[1]:
            new_tiff_data.append(single_data)
    meta['OME']['Image']['Pixels']['TiffData'] = new_tiff_data

    new_config = "4"
    for tile in meta['OME']['Image']['ca:CustomAttributes']['TileConfiguration']['@TileConfiguration'].split("  ")[1:]:
        print(tile.split("Z")[1][:4])
        slice_num = int(tile.split("Z")[1][:4])/slice_step
        if(slice_num < za.darr.shape[1]):
            new_config += "  " + tile
    meta['OME']['Image']['ca:CustomAttributes']['TileConfiguration']['@TileConfiguration'] = new_config
    return meta


def output_downsampled_tiff(output, members_list, za, meta, slice_step):
    """
    Create the new tiff files with the downsampled data and updated 
    metadata
    """
    with tarfile.open(output, 'w') as tar:
        for member in members_list:
            member_slice = int(int(member.name.split("Z")[1][:4])/slice_step)
            channel = int(member.name.split("C")[1][:2])
            if(member_slice < za.darr.shape[1]):
                if(member_slice == 0):
                    new_description = xmltodict.unparse(meta)
                    new_description = new_description.encode("UTF-8")
                    with tifffile.TiffWriter(member.name) as tw:
                        new_data = np.array(za.darr)[channel, member_slice,:,:]
                        tw.write(new_data, description=new_description, metadata=None, planarconfig="CONTIG")
                else:
                    with tifffile.TiffWriter(member.name) as tw:
                        new_data = np.array(za.darr)[channel, member_slice,:,:]
                        tw.write(new_data, metadata=None, planarconfig="CONTIG")
                tar.add(member.name, arcname=member.name)
            os.remove(member.name)

@app.command()
def complete_tiff_downsampling(path_to_source_tar:Annotated[str, typer.Argument(help="ex: dir1/dir2/dataset.tar")],
                                path_to_output_tar:Annotated[str, typer.Argument(help="ex: dir1/dir2/test_dataset.tar")],
                                  ds_x: int=1, ds_y: int=1, ds_z: int=1, slice_step: int=1):
    """
    Make executable from command line using typer commands
    """
    meta, data, member_list = downsample_tiff(path_to_source_tar, ds_x, ds_y, ds_z, slice_step)
    meta = basic_meta_update(meta, data, ds_x,ds_y,ds_z)
    meta = advanced_meta(meta, data, slice_step)
    output_downsampled_tiff(path_to_output_tar, member_list, data, meta, slice_step)
    return meta



if __name__ == "__main__":
    app()