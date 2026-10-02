"""Demonstrates two gotchas with backed=('r') AnnData views (common on
these files, which run gigabytes each) and the working pattern for each.

Usage:
    python backed_anndata_gotchas.py path/to/file.h5ad
"""
import sys

import numpy as np
import anndata as ad


def demo(h5ad_path: str):
    a = ad.read_h5ad(h5ad_path, backed="r")
    view = a[:50]  # a backed *view*, nothing materialized yet

    # Gotcha 1: a backed view can't be .copy()'d directly.
    try:
        view.copy()
        print("(unexpected: .copy() succeeded on this anndata version)")
    except Exception as e:
        print(f"view.copy() raised: {type(e).__name__}: {e}")
        print("  -> anndata's suggested fix is .to_memory(), but that eagerly pulls")
        print("     ALL of obs/uns/varm/obsp too, which can be slow if the file has")
        print("     many small per-key arrays in .uns (these files often do).")

    # Working pattern: materialize only the fields you actually need.
    X = view.X
    X = np.asarray(X.todense() if hasattr(X, "todense") else X[...], dtype=np.float32)
    coords = np.asarray(view.obsm["spatial_3D"]) if "spatial_3D" in view.obsm else None
    print(f"materialized only what's needed: X{X.shape}"
          + (f", spatial_3D{coords.shape}" if coords is not None else ", no spatial_3D")
          + ", obs/uns/varm left untouched")

    # Gotcha 2: row-slicing a CSC-stored matrix. A CSC (compressed-sparse-column)
    # matrix stores data column-major, so any row access touches every column's
    # data structure, and scattered row indices are slow. A contiguous block is
    # faster but is one spatial region of the embryo (README, gotcha 7), so do
    # not use it for anything you submit or score.
    print(f"\nX format on disk: {type(a.X).__name__ if not hasattr(a.X, 'format') else a.X.format}")
    print("  if this says 'csc': sample rows at random, sort the indices first")
    print("  (a[np.sort(idx)]) or read X into memory once and sample from that.")
    print("  Avoid a[start:start+n] for submissions: rows are in spatial order.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    demo(sys.argv[1])
