import os
import tarfile
import tifffile
import xmltodict
import typer
from typing_extensions import Annotated
app = typer.Typer()

def get_members_hybrid(source_dir, step, start_tile_x, start_tile_y, num_tiles_x, num_tiles_y):
    """
    Opens up the given dataset in the tarfile and will extract only the needed members
    based ont he slice step and tiles chosen by the user.
    """
    members_wanted = []
    pairs_wanted = []
    slices_wanted = []

    # get list of all wanted pairs
    for numx in range(start_tile_x, start_tile_x+num_tiles_x):
        numx = "0"+str(numx)
        for numy in range(start_tile_y, start_tile_y+num_tiles_y):
            numy = "0"+str(numy)
            pair = (numx,numy)
            pairs_wanted.append(pair)
    
    with tarfile.open(source_dir, 'r') as tar:
        members = tar.getmembers()
        slice = 0
        
        # add potential slice numbers to a list
        # will go quite over and could use fixing 
        while(slice < len(members)):
            slices_wanted.append(slice)
            slice+=step
        max_slice = 0
        for member in members:
            # get the pair and the slice of each member in the tarfile
            grid = member.name.split("[")[1][:7]
            x_num =grid.split(" x ")[0]
            y_num = grid.split(" x ")[1]
            pair = (x_num,y_num)
            slice = int(member.name.split("Z")[1][:4])

            # find max slice
            if(slice>max_slice):
                max_slice = slice

            # extract the correct tiff files from tar file
            if(pair in pairs_wanted and slice in slices_wanted):
                members_wanted.append(member)
                tar.extract(member)

    # Remove slices that are greater than the max slice
    slices_wanted = (slice for slice in slices_wanted if slice<= max_slice)
    slices_wanted = list(slices_wanted)
    return members_wanted, slices_wanted, pairs_wanted


def correct_metadata_tile(members, slices_wanted, pairs_wanted):
    """
    Goes through each member in the given tar file and if it is the 0th slice 
    which contains all the data and metadata. It will adjust it to the configuartions set, such as
    what slices should be in th enew dataset as well as what tiles to include
    """
    for member in members:
        # get the slice and the channel for each member
        member_slice = int(member.name.split("Z")[1][:4])
        channel = int(member.name.split("C")[1][:2])
        if(member_slice == 0):
            # load tiff file and access metadata
            with tifffile.TiffFile(member.name) as tif:
                # update data to only include the files own data
                data = tif.series[0].asarray()
                data = data[channel][member_slice][:][:]
                meta_dict = xmltodict.parse(tif.ome_metadata)

                # Set to the correct number of slices
                meta_dict['OME']['Image']['Pixels']['@SizeZ'] = len(list(slices_wanted))

                # Remove excess tile configurations from metadata based on tiles
                tile_files = meta_dict['OME']['Image']['ca:CustomAttributes']['TileConfiguration']['@TileConfiguration'].split("  ")[1:]
                proper_files = []
                for file in tile_files:
                    grid = file.split("[")[1][:7]
                    x_num = grid.split(" x ")[0]
                    y_num = grid.split(" x ")[1]
                    pair = (x_num, y_num)
                    if(pair in pairs_wanted):
                        proper_files.append(file)
                new_tile_config = "4"        
                for file in proper_files:
                    new_tile_config += "  " + file
                meta_dict['OME']['Image']['ca:CustomAttributes']['TileConfiguration']['@TileConfiguration'] = new_tile_config
                
                # Remove excess tile configurations from metadata based on slice
                tile_files = meta_dict['OME']['Image']['ca:CustomAttributes']['TileConfiguration']['@TileConfiguration'].split("  ")[1:]
                proper_files = []
                for file in tile_files:
                    if(int(file.split("Z")[1][:4]) in slices_wanted):
                        proper_files.append(file)
                new_tile_config = "4"        
                for file in proper_files:
                    new_tile_config += "  " + file
                meta_dict['OME']['Image']['ca:CustomAttributes']['TileConfiguration']['@TileConfiguration'] = new_tile_config
                
                # Remove excess TiffData from metadata
                tiff_data = meta_dict['OME']['Image']['Pixels']['TiffData']
                new_tiff_data = []
                for single_data in tiff_data:
                    filename = single_data['UUID']['@FileName']
                    grid = filename.split("[")[1][:7]
                    x_num = grid.split(" x ")[0]
                    y_num = grid.split(" x ")[1]
                    pair = (x_num, y_num)
                    if(pair in pairs_wanted):
                        new_tiff_data.append(single_data)
                meta_dict['OME']['Image']['Pixels']['TiffData'] = new_tiff_data
                tiff_data = meta_dict['OME']['Image']['Pixels']['TiffData']
                new_tiff_data = []
                for single_data in tiff_data:
                    if(int(single_data['@FirstZ']) in list(slices_wanted)):
                        new_tiff_data.append(single_data)
                meta_dict['OME']['Image']['Pixels']['TiffData'] = new_tiff_data
                
                # write out new metadata and adjusted data
                new_description = xmltodict.unparse(meta_dict)
                new_description = new_description.encode("UTF-8")
            with tifffile.TiffWriter(member.name) as tw:
                tw.write(data, description=new_description, metadata=None, planarconfig="CONTIG")

def make_tar(members, output):
    """
    Opens up a new tarfile at the specified output and will add all the members that
    were selected from the main dataset and remove the local tiff files produced
    """
    with tarfile.open(output, 'w') as tar:
        for member in members:
            tar.add(member.name, arcname=member.name)
            os.remove(member.name)

@app.command()
def create_test_subset_hybrid(path_to_source_tar:Annotated[str, typer.Argument(help="ex: dir1/dir2/dataset.tar")],
                              path_to_output_tar:Annotated[str, typer.Argument(help="ex: dir1/dir2/tes_dataset.tar")],
                              slice_step: int=20, x_start: int=0,num_x: int=3,y_start: int=0,num_y: int=3):
    members, wanted_slices, wanted_pairs = get_members_hybrid(path_to_source_tar,slice_step, x_start,y_start, num_x, num_y)
    correct_metadata_tile(members, wanted_slices, wanted_pairs)
    make_tar(members, path_to_output_tar)


if __name__ == "__main__":
    app()