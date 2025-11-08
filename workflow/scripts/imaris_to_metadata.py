import json
import re
import h5py
import xmltodict

# -----------------------------
# Helper functions
# -----------------------------

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
        "NumericalAperture": float(custom_attrs.get("ObjectiveNA", {}).get("@ObjectiveNA", 0.0)),
        "Magnification": magnification,
        "OtherAcquisitionParameters": custom_attrs.get("MeasurementMode", {}).get("@MeasurementMode"),
        "InstrumentModel": custom_attrs.get("InstrumentMode", {}).get("@InstrumentMode"),
        "SoftwareVersions": custom_attrs.get("ImspectorVersion", {}).get("@ImspectorVersion"),
    }

    # --- Collect non-BIDS fields into ExtraMetadata ---
    excluded = {
        "ObjectiveMedium", "ObjectiveNA", "ObjectiveID", "MeasurementMode",
        "InstrumentMode", "ImspectorVersion", "DataAxis0", "DataAxis1", "DataAxis2", "AxesLabels"
    }
    extra_metadata = {k: clean_keys(v) for k, v in custom_attrs.items() if k not in excluded}
    bids_json["ExtraMetadata"] = extra_metadata

    return bids_json


# -----------------------------
# Main extraction
# -----------------------------

with h5py.File(snakemake.input.ims, "r") as hdf5_file:
    xml_data = hdf5_file["DataSetInfo/OME Image Tags/Image 0"][:]

xml_dict = parse_xml_bytes(xml_data)
if not xml_dict:
    raise ValueError("Failed to parse XML from .ims file")

custom_attrs = xml_dict["root"]["ca:CustomAttributes"]
bids_metadata = build_bids_metadata(custom_attrs)

# -----------------------------
# Write to JSON
# -----------------------------
with open(snakemake.output.metadata_json, "w", encoding="utf-8") as fp:
    json.dump(bids_metadata, fp, indent=4, ensure_ascii=False)


