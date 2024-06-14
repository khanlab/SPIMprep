import os
import json
import math
from distutils.dir_util import copy_tree
from ome_zarr.io import parse_url
from ome_zarr.reader import Reader
import matplotlib.pyplot as plt
import numpy as np

# arguments for creating the flatfield correction comparisons
ff_s_start=snakemake.params.ff_s_start
ff_s_step=snakemake.params.ff_s_step
ff_cmap=snakemake.params.ff_cmap
ff_corr = snakemake.input.corr
ff_uncorr = snakemake.input.uncorr

# arguments for creating whole slice images
ws_s_step=snakemake.params.ws_s_step
ws_s_start=snakemake.params.ws_s_start
ws_cmap=snakemake.params.ws_cmap
ome_zarr= snakemake.input.ome

# Get output files
output = snakemake.output.out
out_dir = output.split("/")[1]

def make_directories():
   '''
   This function produces a directory for the subject to be able to view
   the proper images 
   '''
   try:
      os.mkdir(f"qc_viewer/{out_dir}")
   except:
      pass
   
   #copy viewers into their respective directory
   copy_tree("qc_viewer/resources/volumeViewer", f"qc_viewer/{out_dir}/volumeViewer")
   copy_tree("qc_viewer/resources/sliceViewer", f"qc_viewer/{out_dir}/sliceViewer")   

def produce_ff_images(corrected_path, uncorrected_path, colour="gray", slice_start=0, slice_step=1):
    '''
    This function produces an html file containing images for each tile before and after flatfield correction
    It can be decided how many slices are wanted.
    '''

    # Get corrected and uncorrected file paths
    corr_zarr = corrected_path
    uncorr_zarr = uncorrected_path

    # Read the corrected and uncorrected ome-zarr data
    proc_reader= Reader(parse_url(corr_zarr))
    unproc_reader = Reader(parse_url(uncorr_zarr))

    proc_data=list(proc_reader())[0].data
    unproc_data = list(unproc_reader())[0].data

    # create a root folder to store the corrected and uncorrected images
    root = "root"
    make_directories()
    os.makedirs(f"qc_viewer/{out_dir}/sliceViewer/{root}/images/corr", exist_ok=True)
    os.mkdir(f"qc_viewer/{out_dir}/sliceViewer/{root}/images/uncorr")

    # Create html file for flatfield corrected images
    with open(f"qc_viewer/{out_dir}/sliceViewer/ff_corr.html", 'w') as fw:
        fw.write(f'''
<!DOCTYPE html>
<head>
  <title>FlatField Correction Check</title>
  <link rel="stylesheet" href="style.css"/>
</head>
<body>
  <a href="../index.html">Back</a>
  <label class="expand-options">
    Expand Images on Click
    <input type="checkbox" id="expand">
    <input type="number" id="expand_scale" value="2">
  </label>
  <table id="table">
    <tbody>
      <h1>Before and After Flatfield Correction</h1>                 
''')
        
        # Add images into the table for each chunk and channel
        for chunk,(tile_corr, tile_uncorr) in enumerate(zip(proc_data[0], unproc_data[0])):
            for chan_num, (channel_corr, channel_uncorr) in enumerate(zip(tile_corr, tile_uncorr)):
                slice = slice_start
                fw.write(f"""      <tr>
        <td colspan={len(channel_corr/slice_step)*2}>
          <h2>Chunk - {chunk}  Channel - {chan_num}</h2>
        </td>
      </tr>
      <tr>""")
                # Process every wanted slice within a chunk and channe;
                while(slice<len(channel_corr)):
                    # Get the contrast limits by removing the absolute highest and lowest data
                    sorted_array = np.sort(np.array(channel_corr[slice].flatten()))[::-1]
                    cmax = sorted_array[math.floor(len(sorted_array)*1/100)]
                    cmin = sorted_array[math.floor(len(sorted_array)*99/100)]
                    
                    # clip and plot data then save image to wanted file path
                    clipped_data_corr = np.clip(channel_corr[slice],cmin,cmax)
                    clipped_data_uncorr = np.clip(channel_uncorr[slice], cmin, cmax)
                    corrected_img_path = f"{root}/images/corr/chunk-{chunk}-channel-{chan_num}-slice-{slice}.jpg"
                    uncorrected_img_path = f"{root}/images/uncorr/chunk-{chunk}-channel-{chan_num}-slice-{slice}.jpg"
                    plt.imsave(f"qc_viewer/{out_dir}/sliceViewer/"+corrected_img_path, clipped_data_corr, cmap=colour)
                    plt.imsave(f"qc_viewer/{out_dir}/sliceViewer/"+uncorrected_img_path, clipped_data_uncorr, cmap=colour) 

                    # Add the images into the html format                  
                    fw.write(f'''       
        <td>
            <img src={corrected_img_path}></img>
            <h3>Corrected</h3>
            <p>Slice-{slice}</p>
        </td>
        <td>
            <img src={uncorrected_img_path}></img>
            <h3>Uncorrected</h3>
            <p>Slice-{slice}</p>
        </td>''')
                    # Increase by user given slice step
                    slice += slice_step

                fw.write("      </tr>")

        fw.write("""
      </tr>
    </tbody>
  </table>
  <script src=image_expand.js></script>                            
                 """)
                

def produce_whole_slice_images(proc_path, colour="gray", slice_start=0, slice_step=1):
    '''
    This function produces an html page containing images of each slice for all channels.
    The slice step can be chosen in the config to choose how many slices you want to view
    '''

    # read ome-zarr data and convert to list
    proc_reader= Reader(parse_url(proc_path))
    proc_data=list(proc_reader())[0].data

    # dump the list into a json to be volume rendered
    produce_json(proc_data)

    # create another root directory to hold the fully preprocessed slice images
    root = "root"
    os.makedirs(f"qc_viewer/{out_dir}/sliceViewer/{root}/images", exist_ok=True)

    # create html page for the whole slices 
    with open(f"qc_viewer/{out_dir}/sliceViewer/whole_slices.html", 'w') as fw:
        fw.write('''
<!DOCTYPE html>
<head>
  <title>Processed Slices</title>
  <link rel="stylesheet" href="style.css"/>
</head>
<body>
  <a href="../index.html">Back</a>
  <label class="expand-options">
    Expand Images on Click
    <input type="checkbox" id="expand">
    <input type="number" id="expand_scale" value="2">
  </label>
  <table id="table">
    <tbody>
      <h1>Processed Image Slices</h1>                 
''')
        # for each channel add the images of the most downsampled data
        for chan_num, channel in enumerate(proc_data[-1]):
            slice = slice_start
            fw.write(f"""      <tr>
      <td colspan={len(channel/slice_step)*2}>
        <h2>Channel - {chan_num}</h2>
      </td>
    </tr>
    <tr>""")
            num_images = 0
            while(slice<len(channel)):
                # Get contrast limits
                sorted_array = np.sort(np.array(channel[slice].flatten()))[::-1]
                cmax = sorted_array[math.floor(len(sorted_array)*1/100)]
                cmin = sorted_array[math.floor(len(sorted_array)*99/100)]

                # normalize data and save image 
                normalize_data= np.clip(channel[slice],cmin,cmax)
                proc_img_path = f"{root}/images/channel-{chan_num}-slice-{slice}.jpg"
                plt.imsave(f"qc_viewer/{out_dir}/sliceViewer/"+proc_img_path, normalize_data, cmap=colour)

                # for every 5th image create a new row in the table
                if(num_images%5 == 0):
                    new_html_text = f"""
      </tr>
      <tr>
        <td>
          <img src="{proc_img_path}"></img>
          <p>Slice-{slice}</p>
        </td>                   
"""
                else:
                    new_html_text = f"""
        <td>
          <img src="{proc_img_path}"></img>
          <p>Slice-{slice}</p>
        </td>
"""
                fw.write(new_html_text)
                slice += slice_step
                num_images+=1

        # end html table
        fw.write("""
      </tr>
    </tbody>
  </table>
  <script src="image_expand.js"></script>   
  </body>                         
                 """)

def produce_json(data):
    '''
    Produces a json file containing the most downsampled
    image data to be volume rendered into a 3D image.
    '''

    with open(f"qc_viewer/{out_dir}/volumeViewer/volumeData.json", 'w') as f:
        data = np.array(data[-1]).tolist()
        data = json.dumps(data)
        f.write(data)


def combine_sample_htmls(ffcorr_html, proc_html):
    '''
    Produces and index.html page connecting the two image reports as well as 
    the 3D volume rendering page
    '''

    with open(f'qc_viewer/{out_dir}/index.html', 'w') as f:
        f.write(f"""
<!DOCTYPE html>
<head>
  <title>Processed Slices</title>
  <link rel="stylesheet" href="style.css"/>
</head>
<body>
  <a href="../index.html">Back</a>
  <h1>{out_dir.split("-")[0]}</h1>
  <a href="sliceViewer/{ffcorr_html}">Flatfield Correction Before and After</a>
  <br>
  <a href="sliceViewer/{proc_html}">Full Processed Slices</a>
  <br>
  <a href="volumeViewer/volRender.html">3D Image</a>
<body>       
                """)        


def create_main_html():
   """
   This function creates an html file connecting all the samples viewers together.
   If the file is empty it will produce the header and if it is not it just adds 
   another sampel link
   """

   file="qc_viewer/index.html"
   
   with open(file, 'a') as f:
      if(os.path.getsize(file) <= 20):
         f.write(f"""
<!DOCTYPE html>
<head>
  <title>Sample Check</title>
  <link rel="stylesheet" href="style.css"/>
</head>
<body>
  <h1>Subject Reports</h1>
  <a href="{out_dir}/index.html">{out_dir.split('-')[0]}</a>
  <br>
                 """)
      else:
         f.write(f"""
  <a href="{out_dir}/index.html">{out_dir.split('-')[0]}</a>
  <br>
                 """)
      


produce_ff_images(ff_corr, ff_uncorr, slice_start=ff_s_start, slice_step=ff_s_step,
                  colour=ff_cmap)
produce_whole_slice_images(ome_zarr, slice_start=ws_s_start, slice_step=ws_s_step,
                            colour=ws_cmap)
combine_sample_htmls("ff_corr.html", "whole_slices.html")
create_main_html()



