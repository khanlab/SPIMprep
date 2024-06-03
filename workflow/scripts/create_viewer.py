import matplotlib.pyplot as plt
import numpy as np
from ome_zarr.io import parse_url
from ome_zarr.reader import Reader
import os
import json
import math
from distutils.dir_util import copy_tree

generate_report = snakemake.params.create_report

ff_s_start=snakemake.params.ff_slice_start
ff_s_step=snakemake.params.ff_slice_step
ff_cmap=snakemake.params.ff_cmap
ff_corr = snakemake.input.corr
ff_uncorr = snakemake.input.uncorr

ws_s_step=snakemake.params.ws_s_step
ws_s_start=snakemake.params.ws_s_start
ws_cmap=snakemake.params.ws_cmap

ome_zarr= snakemake.input.ome

output = snakemake.output.out
out_dir = output.split("/")[1]

def make_directories():
   '''
   This function produces a directory for the sample to be able to view
   the proper images 
   '''
   try:
      os.mkdir(f"qc_viewer/{out_dir}")
   except:
      pass
   
   #copy viewers into the directories
   copy_tree("qc_viewer/resources/volumeViewer", f"qc_viewer/{out_dir}/volumeViewer")
   copy_tree("qc_viewer/resources/sliceViewer", f"qc_viewer/{out_dir}/sliceViewer")   

def produce_ff_images(corrected_path, uncorrected_path, colour="gray", slice_start=0, slice_step=1):
    '''
    This function produces an html file containing images for each tile before and after flatfield correction
    It can be decided how many slices are wanted. The contrast can also be adjusted for 
    different images and channels to ensure good quality.
    '''
    corr_zarr = corrected_path
    uncorr_zarr = uncorrected_path

    proc_reader= Reader(parse_url(corr_zarr))
    unproc_reader = Reader(parse_url(uncorr_zarr))

    proc_data=list(proc_reader())[0].data
    unproc_data = list(unproc_reader())[0].data
    root = "root"
    make_directories()
    os.makedirs(f"qc_viewer/{out_dir}/sliceViewer/{root}/images/corr", exist_ok=True)
    os.mkdir(f"qc_viewer/{out_dir}/sliceViewer/{root}/images/uncorr")
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
        
        for chunk,(tile_corr, tile_uncorr) in enumerate(zip(proc_data[0], unproc_data[0])):
            for chan_num, (channel_corr, channel_uncorr) in enumerate(zip(tile_corr, tile_uncorr)):
                slice = slice_start
                fw.write(f"""      <tr>
        <td colspan={len(channel_corr/slice_step)*2}>
          <h2>Chunk - {chunk}  Channel - {chan_num}</h2>
        </td>
      </tr>
      <tr>""")
                while(slice<len(channel_corr)):
                    """
                    if(chan_num==0):
                        cmin = ch1_contrast_low
                        cmax= ch1_contrast_high
                    else:
                        cmin=ch2_contrast_low
                        cmax=ch2_contrast_high
                    """
                    sorted_array = np.sort(np.array(channel_corr[slice].flatten()))[::-1]
                    cmax = sorted_array[math.floor(len(sorted_array)*1/100)]
                    cmin = sorted_array[math.floor(len(sorted_array)*99/100)]
                    print(f"chan-{chan_num} slice-{slice} min-{cmin} max-{cmax}")
                    #cmin = np.amin(channel_corr[slice])
                    #cmax = np.amax(channel_uncorr[slice])
                    normalize_data_corr = np.clip(channel_corr[slice],cmin,cmax)
                    normalize_data_uncorr = np.clip(channel_uncorr[slice], cmin, cmax)
                    corrected_img_path = f"{root}/images/corr/chunk-{chunk}-channel-{chan_num}-slice-{slice}.jpg"
                    uncorrected_img_path = f"{root}/images/uncorr/chunk-{chunk}-channel-{chan_num}-slice-{slice}.jpg"
                    plt.imsave(f"qc_viewer/{out_dir}/sliceViewer/"+corrected_img_path, normalize_data_corr, cmap=colour)
                    plt.imsave(f"qc_viewer/{out_dir}/sliceViewer/"+uncorrected_img_path,normalize_data_uncorr, cmap=colour)                   
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
    This function produces an html page containing images of each slice for both channels.
    The slice step can be chosen as well as the contrast for each channel to ensure the best quality images
    '''
    proc_reader= Reader(parse_url(proc_path))
    proc_data=list(proc_reader())[0].data
    produce_json(proc_data)

    root = "root"
    os.makedirs(f"qc_viewer/{out_dir}/sliceViewer/{root}/images", exist_ok=True)
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
        for chan_num, channel in enumerate(proc_data[5]):
            slice = slice_start
            fw.write(f"""      <tr>
      <td colspan={len(channel/slice_step)*2}>
        <h2>Channel - {chan_num}</h2>
      </td>
    </tr>
    <tr>""")
            num_images = 0
            while(slice<len(channel)):
                """
                if(chan_num==0):
                    cmin = ch1_contrast_low
                    cmax= ch1_contrast_high
                else:
                    cmin=ch2_contrast_low
                    cmax=ch2_contrast_high
                    """
                sorted_array = np.sort(np.array(channel[slice].flatten()))[::-1]
                cmax = sorted_array[math.floor(len(sorted_array)*1/100)]
                cmin = sorted_array[math.floor(len(sorted_array)*99/100)]
                print(f"{chan_num} {cmax} {cmin}")
                normalize_data= np.clip(channel[slice],cmin,cmax)
                proc_img_path = f"{root}/images/channel-{chan_num}-slice-{slice}.jpg"
                plt.imsave(f"qc_viewer/{out_dir}/sliceViewer/"+proc_img_path, normalize_data, cmap=colour)
                new_html_text = ""
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
        fw.write("""
      </tr>
    </tbody>
  </table>
  <script src="image_expand.js"></script>   
  </body>                         
                 """)

def produce_json(data):
    '''
    Produces a json format containing the most downsampled
    image data to be rendered into a 3D image.
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
  <br>
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
      

if(generate_report):
  """
  Runs entire script if is configured to in the config.yml file
  If not it will just produce the target file for rule to pass
  """
  produce_ff_images(ff_corr, ff_uncorr, slice_start=ff_s_start, slice_step=ff_s_start,
                    colour=ff_cmap)
  produce_whole_slice_images(ome_zarr, slice_start=ws_s_start, slice_step=ws_s_step,
                             colour=ws_cmap)
  combine_sample_htmls("ff_corr.html", "whole_slices.html")
  create_main_html()
else:
   directories=f"qc_viewer/{out_dir}/volumeViewer"
   os.makedirs(directories, exist_ok=True)
   with open(directories+"/volumeData.json", 'w') as f:
      f.write("")


