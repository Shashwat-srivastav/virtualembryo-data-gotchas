"""Check whether an .h5ad file's obs_names (cell barcodes) are actually unique.
At least one released file in this Challenge has duplicated barcodes from a
segmentation-pipeline artifact (same slide/section, different internal
sample_id) -- joining on obs_names there will merge unrelated cells.

Usage:
    python check_duplicate_obs_names.py path/to/file.h5ad
"""
import sys
from collections import Counter

import anndata as ad


def check(h5ad_path: str) -> bool:
    a = ad.read_h5ad(h5ad_path, backed="r")
    names = list(a.obs_names)
    counts = Counter(names)
    dups = {n: c for n, c in counts.items() if c > 1}

    if not dups:
        print(f"OK: {h5ad_path} -- all {len(names)} obs_names unique")
        return True

    n_dup_cells = sum(dups.values())
    print(f"NOT UNIQUE: {h5ad_path} -- {len(dups)} barcodes repeated, "
          f"{n_dup_cells} of {len(names)} cells affected")
    example = next(iter(dups))
    print(f"  example: '{example}' appears {dups[example]} times")
    print("  do not join or index on obs_names for this file -- use positional "
          "indices, or a per-file unique key if one exists in .obs, instead")
    return False


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    ok = check(sys.argv[1])
    raise SystemExit(0 if ok else 1)
