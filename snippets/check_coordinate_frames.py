"""Print the spatial_3D bounding box for one or more .h5ad files, to reveal
whether coordinates are on a shared frame or are per-file local (differing in
both scale and origin). Pooling coordinates across files without re-centering
lets a model shortcut-learn on absolute offset instead of shape.

Usage:
    python check_coordinate_frames.py file1.h5ad file2.h5ad ...
"""
import sys

import numpy as np
import anndata as ad


def bbox(h5ad_path: str):
    a = ad.read_h5ad(h5ad_path, backed="r")
    if "spatial_3D" not in a.obsm:
        print(f"{h5ad_path}: no obsm['spatial_3D']")
        return
    C = np.asarray(a.obsm["spatial_3D"])
    lo, hi = C.min(axis=0), C.max(axis=0)
    centroid = C.mean(axis=0)
    print(f"{h5ad_path}: n={C.shape[0]}  "
          f"x in [{lo[0]:.1f},{hi[0]:.1f}] y in [{lo[1]:.1f},{hi[1]:.1f}] z in [{lo[2]:.1f},{hi[2]:.1f}]  "
          f"centroid=({centroid[0]:.1f},{centroid[1]:.1f},{centroid[2]:.1f})")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    for path in sys.argv[1:]:
        bbox(path)
    print("\nIf these boxes don't share a common scale/origin, re-center (and "
          "optionally rescale) each file to its own centroid before pooling "
          "coordinates across files for any model.")
