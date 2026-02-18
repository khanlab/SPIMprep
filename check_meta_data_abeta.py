import json
from pathlib import Path

base = Path("/nfs/khan/datasets/selma_all/Ab_plaques/bids")

for zarr_path in sorted(base.glob("sub-*/micr/*.ome.zarr")):
    print(f"\n{'='*60}")
    print(f"Zarr: {zarr_path.name}")
    
    # Root zarr.json
    root_meta = zarr_path / "zarr.json"
    if root_meta.exists():
        with open(root_meta) as f:
            meta = json.load(f)
        print(json.dumps(meta, indent=2)[:3000])
    
    # Check each level's zarr.json
    for level_dir in sorted(zarr_path.iterdir()):
        if level_dir.is_dir():
            level_meta = level_dir / "zarr.json"
            if level_meta.exists():
                with open(level_meta) as f:
                    lmeta = json.load(f)
                print(f"\n  Level '{level_dir.name}': shape={lmeta.get('shape')}, chunks={lmeta.get('chunk_grid', {})}, dtype={lmeta.get('data_type')}")
