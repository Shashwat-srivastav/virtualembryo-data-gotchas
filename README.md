# Virtual Embryo Challenge — Data Gotchas

Six things about the Challenge's released `.h5ad` files that cost real debugging
time to find, written up so you don't have to rediscover them. Each one has a
one-paragraph explanation and a runnable check in [`snippets/`](snippets/) you
can point at your own downloaded copy of the data.

This isn't about modeling strategy — it's data engineering: how the files are
actually stored, and where a plausible-looking naive assumption breaks.

## 1. Gene panels can differ in size across stages

Some early-stage files measure fewer genes than the full panel used elsewhere
(the T2-embryo panel drops two genes that later stages do have). If you assume
every file uses the same gene set and just index by position, you'll silently
misalign columns instead of getting an error. Always verify `var_names`
against the official `*.genes.txt` panel for that board, byte-for-byte and in
order, before treating a file's columns as aligned to anything else. And note
the corollary: a prediction can't be scored against a gene the *reference*
stage never observed, even if the *target* stage did measure it.

```bash
python snippets/check_gene_panel.py path/to/file.h5ad path/to/panel.genes.txt
```

## 2. `obs_names` are not guaranteed unique

At least one released file has barcodes that appear twice — same slide, same
section, different internal `sample_id` — a segmentation-pipeline artifact,
not a duplicate measurement. If you join or index anything on `obs_names`
assuming uniqueness, you will silently merge two different cells under one
key. Check first; use positional indices (or another per-file unique key) if
you find duplicates.

```bash
python snippets/check_duplicate_obs_names.py path/to/file.h5ad
```

## 3. Spatial coordinates are per-file local, not globally registered

`spatial_3D` isn't just "unregistered" between files — it differs in **scale
and origin**, not only rotation. Some files sit near the origin with a small
bounding box; others sit at a completely different offset with a much larger
one. If you pool coordinates across files without re-centering each one
first, a model can shortcut-learn on absolute position instead of the actual
spatial structure you want it to learn.

```bash
python snippets/check_coordinate_frames.py file1.h5ad file2.h5ad ...
```

## 4. Backed AnnData views can't be `.copy()`'d directly

Opening a large file with `backed="r"` (sensible — these files run gigabytes)
gives you a view, and calling `.copy()` on that view raises an error asking
you to use `.to_memory()` instead. The catch: `.to_memory()` eagerly pulls
**everything** — all of `.obs`, `.uns`, `.varm` — not just the slice you
wanted, which can be slow when a file's `.uns` holds many small per-key
arrays (common in these files). If you only need `X` and maybe `obsm`,
materialize those fields directly instead of calling `.to_memory()` on the
whole object.

```bash
python snippets/backed_anndata_gotchas.py path/to/file.h5ad
```

## 5. `X` is stored CSC, not dense

Row (cell) slicing on a CSC (compressed-sparse-**column**) matrix touches
every column's data structure no matter what you do — but a **contiguous**
row block is still far cheaper to materialize than a **scattered** random
subsample of the same size, because each column's slice can stay a single
contiguous range instead of becoming a fancy-index gather. If your use case
doesn't specifically need a random subsample (e.g. you're building a floor
baseline that resubmits an observed stage "verbatim" anyway), prefer
`view = a[start:start+n]` over `view = a[random_indices]`.

Covered together with gotcha 4 in the same script:

```bash
python snippets/backed_anndata_gotchas.py path/to/file.h5ad
```

## 6. Cell-type labels aren't harmonized across stages

The same biological population can appear under different label spellings at
different stages (for example, one label at an earlier stage split into
several numbered variants later). If a `groupby` or set-diff on cell type
comes back with a label "missing" at some stage, check whether it's actually
absent biologically or just renamed — don't assume the former.

```python
import anndata as ad
a1 = ad.read_h5ad("stage1.h5ad", backed="r")
a2 = ad.read_h5ad("stage2.h5ad", backed="r")
only_in_1 = set(a1.obs["celltype"]) - set(a2.obs["celltype"])
only_in_2 = set(a2.obs["celltype"]) - set(a1.obs["celltype"])
print("only in stage1:", only_in_1)
print("only in stage2:", only_in_2)
# eyeball these for renames (e.g. "X" vs "X-1"/"X-2") before treating either
# side as a true biological absence
```

## Bonus: unaudited precomputed embeddings

Some files ship precomputed embeddings (UMAP/PCA) alongside raw `X`. Their
exact provenance relative to any external-data rules in the Challenge isn't
documented, so treat them as untrusted and prefer recomputing from raw `.X`
yourself rather than consuming them directly.

## Once you're past these: a starter baseline

If you want a working, submission-ready starting point after handling the
gotchas above, see the companion contribution:
**virtualembryo-baseline-celltype-interp** <!-- TODO: link once published -->
— a cell-type-conditioned interpolation baseline, stronger than the naive
`copy_nearest`/`linear` starting points, with the same submission-schema
validation these snippets assume.

## License

MIT — see [LICENSE](LICENSE).

---

Filed as a Community Contribution to the NeurIPS 2026 Virtual Embryo Challenge.
