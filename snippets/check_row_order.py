"""Are a file's rows stored in spatial order? If so, head(N) or any contiguous row block is a spatial slab.

For each file with obsm["spatial_3D"], prints:
  - the strongest |Spearman correlation| between row index and a principal axis of the coordinates
    (a random order gives about 0; a slab-wise order gives a large value);
  - voxel coverage of the first N rows, and of the worst and median contiguous N-row block, as a share of what a
    random sample of N rows covers (16^3 grid over the file's extent; a representative sample is near 100 %).

Only the coordinates are read, so this is fast even on the large stages.

    python snippets/check_row_order.py data/raw/*.h5ad --n 5000
"""
import argparse

import anndata as ad
import numpy as np
from scipy.stats import spearmanr


def coverage(C: np.ndarray, lo: np.ndarray, hi: np.ndarray, bins: int = 16) -> int:
    v = np.clip(((C - lo) / (hi - lo + 1e-12) * bins).astype(int), 0, bins - 1)
    return len(np.unique(v[:, 0] * bins * bins + v[:, 1] * bins + v[:, 2]))


def check(path: str, n: int, seed: int = 0) -> None:
    a = ad.read_h5ad(path, backed="r")
    if "spatial_3D" not in a.obsm:
        print(f"{path}: no obsm['spatial_3D'], skipped")
        return
    C = np.asarray(a.obsm["spatial_3D"], dtype=np.float64)
    N = len(C)
    n = min(n, N)
    Z = (C - C.mean(0)) @ np.linalg.svd(C - C.mean(0), full_matrices=False)[2].T
    rho = max(abs(spearmanr(np.arange(N), Z[:, k])[0]) for k in range(Z.shape[1]))
    lo, hi = C.min(0), C.max(0)
    rng = np.random.default_rng(seed)
    rand = np.mean([coverage(C[rng.choice(N, n, replace=False)], lo, hi) for _ in range(5)])
    starts = np.linspace(0, N - n, num=min(20, N - n + 1)).astype(int)
    blocks = np.array([coverage(C[s:s + n], lo, hi) for s in starts]) / rand
    print(f"{path}: {N} cells | max |rho(row, axis)| {rho:.2f} | coverage vs random {n}-cell sample: "
          f"first {n} rows {coverage(C[:n], lo, hi) / rand:.0%}, contiguous blocks worst {blocks.min():.0%} "
          f"median {np.median(blocks):.0%}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="+")
    p.add_argument("--n", type=int, default=5000, help="cells per sample (default 5000)")
    args = p.parse_args()
    for f in args.files:
        check(f, args.n)


if __name__ == "__main__":
    main()
