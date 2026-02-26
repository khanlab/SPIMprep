import html
import json
import re

import h5py
import xmltodict

# -----------------------------
# Helper functions
# -----------------------------


def parse_xml_str(xml_str):
    """Convert XML string to dict, decoding HTML entities (e.g. &#181; → µ)."""

    def decode_entities(obj):
        """Recursively decode HTML entities in all strings."""
        if isinstance(obj, dict):
            return {k: decode_entities(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [decode_entities(i) for i in obj]
        elif isinstance(obj, str):
            return html.unescape(obj)
        else:
            return obj

    try:
        parsed = xmltodict.parse(xml_str, namespace_separator=":")
        return decode_entities(parsed)

    except Exception as e:
        print(f"Error parsing XML: {e}")
        return None


def parse_xml_bytes(xml_bytes):
    """Convert raw byte array to parsed XML dict with namespace separator."""
    xml_str = bytes(xml_bytes).decode("utf-8", errors="ignore")
    try:
        return xmltodict.parse(f"<root>{xml_str}</root>", namespace_separator=":")
    except Exception as e:
        print(f"Error parsing XML: {e}")
        return None


def clean_keys(obj):
    """Recursively remove '@' from dict keys."""
    if isinstance(obj, dict):
        return {k.replace("@", ""): clean_keys(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [clean_keys(i) for i in obj]
    return obj


def h5_attr_to_str(attr):
    """Convert an HDF5 byte-array attribute to a Python string.

    Imaris stores attributes as arrays of single-byte values, e.g.
    [b'4', b'.', b'0'] → '4.0'
    """
    if attr is None:
        return None
    return "".join(b.decode("utf-8") for b in attr)


def build_bids_metadata(custom_attrs):
    """Convert custom attribute dict into BIDS microscopy-compliant metadata."""
    # --- PixelSize and Units ---
    px_x = float(custom_attrs["DataAxis0"]["@PhysicalUnit"])
    px_y = float(custom_attrs["DataAxis1"]["@PhysicalUnit"])
    px_z = abs(float(custom_attrs["DataAxis2"]["@PhysicalUnit"]))
    pixel_size = [px_x, px_y, px_z]

    unit_label = custom_attrs.get("AxesLabels", {}).get("@FirstAxis-Unit", "µm")
    unit_label = unit_label.replace("µ", "u") if "µ" in unit_label else unit_label

    # --- Extract Magnification (if present) ---
    obj_id = custom_attrs.get("ObjectiveID", {}).get("@ObjectiveID", "")
    magnification_match = re.search(r"(\d+(?:\.\d+)?)x", obj_id)
    magnification = float(magnification_match.group(1)) if magnification_match else None

    # --- Build BIDS JSON dict ---
    bids_json = {
        "PixelSize": pixel_size,
        "PixelSizeUnits": unit_label,
        "Immersion": custom_attrs.get("ObjectiveMedium", {}).get("@ObjectiveMedium"),
        "NumericalAperture": float(
            custom_attrs.get("ObjectiveNA", {}).get("@ObjectiveNA", 0.0)
        ),
        "Magnification": magnification,
        "OtherAcquisitionParameters": custom_attrs.get("MeasurementMode", {}).get(
            "@MeasurementMode"
        ),
        "InstrumentModel": custom_attrs.get("InstrumentMode", {}).get(
            "@InstrumentMode"
        ),
        "SoftwareVersions": custom_attrs.get("ImspectorVersion", {}).get(
            "@ImspectorVersion"
        ),
    }

    # --- Collect non-BIDS fields into ExtraMetadata ---
    excluded = {
        "ObjectiveMedium",
        "ObjectiveNA",
        "ObjectiveID",
        "MeasurementMode",
        "InstrumentMode",
        "ImspectorVersion",
        "DataAxis0",
        "DataAxis1",
        "DataAxis2",
        "AxesLabels",
    }
    extra_metadata = {
        k: clean_keys(v) for k, v in custom_attrs.items() if k not in excluded
    }
    bids_json["ExtraMetadata"] = extra_metadata

    return bids_json


def build_bids_metadata_from_native_imaris(hdf5_file):
    """Build BIDS metadata directly from native Imaris HDF5 attributes.

    Fallback for .ims files that lack embedded OME XML tags.
    Computes voxel sizes from image extents: pixel_size = (ExtMax - ExtMin) / N_voxels
    """
    img = hdf5_file["DataSetInfo/Image"]

    # --- Image dimensions ---
    nx = int(h5_attr_to_str(img.attrs["X"]))
    ny = int(h5_attr_to_str(img.attrs["Y"]))
    nz = int(h5_attr_to_str(img.attrs["Z"]))

    # --- Physical extents (in file units, typically µm) ---
    ext_min_x = float(h5_attr_to_str(img.attrs["ExtMin0"]))
    ext_max_x = float(h5_attr_to_str(img.attrs["ExtMax0"]))
    ext_min_y = float(h5_attr_to_str(img.attrs["ExtMin1"]))
    ext_max_y = float(h5_attr_to_str(img.attrs["ExtMax1"]))
    ext_min_z = float(h5_attr_to_str(img.attrs["ExtMin2"]))
    ext_max_z = float(h5_attr_to_str(img.attrs["ExtMax2"]))

    # --- Compute voxel sizes ---
    px_x = abs(ext_max_x - ext_min_x) / nx
    px_y = abs(ext_max_y - ext_min_y) / ny
    px_z = abs(ext_max_z - ext_min_z) / nz
    pixel_size = [px_x, px_y, px_z]

    # --- Units ---
    unit_raw = h5_attr_to_str(img.attrs.get("Unit", None))
    if unit_raw is None:
        unit_label = "um"
    else:
        unit_label = unit_raw.replace("µ", "u") if "µ" in unit_raw else unit_raw

    # --- Channel metadata ---
    extra_channels = {}
    ch_idx = 0
    while f"DataSetInfo/Channel {ch_idx}" in hdf5_file:
        ch_grp = hdf5_file[f"DataSetInfo/Channel {ch_idx}"]
        ch_info = {}
        for attr_name in ["Name", "LSMExcitationWavelength", "LSMEmissionWavelength",
                          "Color", "Description"]:
            val = ch_grp.attrs.get(attr_name, None)
            if val is not None:
                ch_info[attr_name] = h5_attr_to_str(val)
        extra_channels[f"Channel_{ch_idx}"] = ch_info
        ch_idx += 1

    # --- Software version ---
    sw_version = None
    if "DataSetInfo/ImarisDataSet" in hdf5_file:
        ds_grp = hdf5_file["DataSetInfo/ImarisDataSet"]
        ver = ds_grp.attrs.get("Version", None)
        if ver is not None:
            sw_version = h5_attr_to_str(ver)

    print(f"Native Imaris metadata: PixelSize={pixel_size}, Unit={unit_label}")
    print(f"  Image dims: X={nx}, Y={ny}, Z={nz}")
    print(f"  Extents X: [{ext_min_x}, {ext_max_x}]")
    print(f"  Extents Y: [{ext_min_y}, {ext_max_y}]")
    print(f"  Extents Z: [{ext_min_z}, {ext_max_z}]")

    bids_json = {
        "PixelSize": pixel_size,
        "PixelSizeUnits": unit_label,
        "Immersion": None,
        "NumericalAperture": 0.0,
        "Magnification": None,
        "OtherAcquisitionParameters": None,
        "InstrumentModel": None,
        "SoftwareVersions": sw_version,
        "ExtraMetadata": {
            "Channels": extra_channels,
            "ImageExtents": {
                "ExtMin": [ext_min_x, ext_min_y, ext_min_z],
                "ExtMax": [ext_max_x, ext_max_y, ext_max_z],
            },
        },
    }

    return bids_json


# -----------------------------
# Main extraction
# -----------------------------


def print_h5_keys_recursive(obj, indent=""):
    """
    Recursively prints the keys and their hierarchy within an h5py object.

    Args:
        obj: An h5py.File or h5py.Group object.
        indent: A string representing the current indentation level for printing.
    """
    if isinstance(obj, h5py.File) or isinstance(obj, h5py.Group):
        for key, value in obj.items():
            print(f"{indent}- {key}")
            if isinstance(value, h5py.Group):
                print_h5_keys_recursive(value, indent + "  ")
    elif isinstance(obj, h5py.Dataset):
        # Datasets are typically leaf nodes, their names are already printed by the parent group
        pass
    else:
        print(f"{indent}- (Unknown HDF5 object type: {type(obj)})")


bids_metadata = None

with h5py.File(snakemake.input.ims, "r") as hdf5_file:

    # ── Strategy 1: OME XML tags embedded in the .ims file ──────────────
    try:
        xml_data = hdf5_file["DataSetInfo/OME Image Tags/Image 0"][:]
        xml_dict = parse_xml_bytes(xml_data)
        if xml_dict:
            custom_attrs = xml_dict["root"]["ca:CustomAttributes"]
            print(custom_attrs)
            bids_metadata = build_bids_metadata(custom_attrs)
    except (KeyError, TypeError, ValueError) as e:
        print(
            f"Warning: cannot find OME metadata from imaris file ({e}), "
            "trying tif fallback..."
        )

    # ── Strategy 2: associated .ome.tif in standard folder structure ────
    if bids_metadata is None:
        try:
            from glob import glob
            from pathlib import Path

            import tifffile

            p = Path(snakemake.input.ims)
            subj = p.parent.name
            raw_root = p.parents[2]
            new_path = str(raw_root / "tif_4x_*" / subj / "*.ome.tif")

            tif_file = sorted(glob(new_path))[0]
            with tifffile.TiffFile(tif_file) as tif:
                xml_data = tif.ome_metadata
                xml_dict = parse_xml_str(xml_data)
                custom_attrs = xml_dict["OME"]["Image"]["ca:CustomAttributes"]
                print(custom_attrs)
                bids_metadata = build_bids_metadata(custom_attrs)
        except (IndexError, KeyError, TypeError, ValueError) as e:
            print(
                f"Warning: cannot find associated tif files ({e}), "
                "trying native Imaris HDF5 attributes..."
            )

    # ── Strategy 3: native Imaris HDF5 attributes (ExtMin/ExtMax) ───────
    if bids_metadata is None:
        try:
            bids_metadata = build_bids_metadata_from_native_imaris(hdf5_file)
        except (KeyError, TypeError, ValueError) as e:
            raise ValueError(
                f"Cannot extract metadata from .ims file using any strategy. "
                f"Native Imaris fallback failed with: {e}"
            ) from e


if bids_metadata is None:
    raise ValueError("Failed to extract metadata from .ims file")


# -----------------------------
# Write to JSON
# -----------------------------
with open(snakemake.output.metadata_json, "w", encoding="utf-8") as fp:
    json.dump(bids_metadata, fp, indent=4, ensure_ascii=False)
