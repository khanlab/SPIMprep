import json

import h5py
import xmltodict

with h5py.File(snakemake.input.ims, "r") as hdf5_file:
    xml_data = hdf5_file["DataSetInfo/OME Image Tags/Image 0"][:]


# Convert byte array to string and then to a dictionary
xml_str = bytes(xml_data).decode(
    "utf-8", errors="ignore"
)  # Decode byte array to string

try:
    xml_dict = xmltodict.parse(f"<root>{xml_str}</root>", namespace_separator=":")
except Exception as e:
    print(f"Error parsing XML: {e}")


metadata = {}
metadata["physical_size_x"] = float(
    xml_dict["root"]["ca:CustomAttributes"]["DataAxis0"]["@PhysicalUnit"]
)
metadata["physical_size_y"] = float(
    xml_dict["root"]["ca:CustomAttributes"]["DataAxis1"]["@PhysicalUnit"]
)
metadata["physical_size_z"] = abs(
    float(xml_dict["root"]["ca:CustomAttributes"]["DataAxis3"]["@PhysicalUnit"])
)
metadata["PixelSize"] = [
    metadata["physical_size_z"] / 1000.0,
    metadata["physical_size_y"] / 1000.0,
    metadata["physical_size_x"] / 1000.0,
]  # zyx since OME-Zarr is ZYX
metadata["PixelSizeUnits"] = "mm"

# write metadata to json
with open(snakemake.output.metadata_json, "w") as fp:
    json.dump(metadata, fp, indent=4)
