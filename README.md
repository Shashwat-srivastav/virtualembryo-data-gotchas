# Virtual Embryo Challenge: Data Gotchas

Seven things about the Challenge's released `.h5ad` files that cost real debugging
time to find, written up so you don't have to rediscover them. Each one has a
one-paragraph explanation and a runnable check in [`snippets/`](snippets/) you
can point at your own downloaded copy of the data. One of them (gotcha 7) cost
us 10.3 leaderboard points before we found it.

Already stuck on an error? Paste it into the bundled
[Gotcha Doctor](#gotcha-doctor-paste-an-error-get-the-gotcha) and it tells you
which of the seven it most likely is.

This is about data engineering, not modeling strategy: how the files are
actually stored, and where a plausible-looking naive assumption breaks.

## 1. Gene panels can differ in size across stages

Some early-stage files measure fewer genes than the full panel used elsewhere
(the T2-embryo panel drops two genes that later stages do have). If you assume
every file uses the same gene set and just index by position, you'll
misalign columns instead of getting an error. Always verify `var_names`
against the official `*.genes.txt` panel for that board, byte-for-byte and in
order, before treating a file's columns as aligned to anything else. And note
the corollary: a prediction can't be scored against a gene the *reference*
stage never observed, even if the *target* stage did measure it.

```bash
python snippets/check_gene_panel.py path/to/file.h5ad path/to/panel.genes.txt
```

## 2. `obs_names` are not guaranteed unique

At least one released file has barcodes that appear twice (same slide, same
section, different internal `sample_id`). This is a segmentation-pipeline artifact,
not a duplicate measurement. If you join or index anything on `obs_names`
assuming uniqueness, you will merge two different cells under one
key. Check first; use positional indices (or another per-file unique key) if
you find duplicates.

```bash
python snippets/check_duplicate_obs_names.py path/to/file.h5ad
```

## 3. Spatial coordinates are per-file local, not globally registered

`spatial_3D` is not only unregistered between files. It differs in **scale
and origin**, not only rotation. Some files sit near the origin with a small
bounding box; others sit at a completely different offset with a much larger
one. If you pool coordinates across files without re-centering each one
first, a model can shortcut-learn on absolute position instead of the actual
spatial structure you want it to learn.

```bash
python snippets/check_coordinate_frames.py file1.h5ad file2.h5ad ...
```

## 4. Backed AnnData views can't be `.copy()`'d directly

Opening a large file with `backed="r"` (sensible, since these files run to gigabytes)
gives you a view, and calling `.copy()` on that view raises an error asking
you to use `.to_memory()` instead. The catch: `.to_memory()` eagerly pulls
**everything** (all of `.obs`, `.uns`, `.varm`), not just the slice you
wanted, which can be slow when a file's `.uns` holds many small per-key
arrays (common in these files). If you only need `X` and maybe `obsm`,
materialize those fields directly instead of calling `.to_memory()` on the
whole object.

```bash
python snippets/backed_anndata_gotchas.py path/to/file.h5ad
```

## 5. `X` is stored CSC, not dense

Row (cell) slicing on a CSC (compressed-sparse-**column**) matrix touches
every column's data structure no matter what you do, and a **contiguous**
row block is far cheaper to materialize than a **scattered** random
subsample of the same size, because each column's slice can stay a single
contiguous range instead of becoming a fancy-index gather.

**Do not let that speed decide which cells you submit or score.** An earlier
version of this page recommended `a[start:start+n]` for floor baselines. That
advice was wrong: rows are stored in spatial order (gotcha 7), so a contiguous
block is one region of the embryo, and it cost us 10.3 leaderboard points on
one board. For anything that will be scored, take a random sample. To keep it
fast, sort the random indices (`np.sort(idx)`) before slicing, or read `X`
into memory once and sample from that.

Covered together with gotcha 4 in the same script:

```bash
python snippets/backed_anndata_gotchas.py path/to/file.h5ad
```

## 6. Cell-type labels aren't harmonized across stages

The same biological population can appear under different label spellings at
different stages (for example, one label at an earlier stage split into
several numbered variants later). If a `groupby` or set-diff on cell type
comes back with a label "missing" at some stage, check whether it's actually
absent biologically or just renamed. Don't assume the former.

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

For candidate renames ranked by string similarity, see the companion
[Cell-Type Vocabulary Diff](https://github.com/Shashwat-srivastav/virtualembryo-celltype-vocab-diff).

## 7. Rows are stored in spatial order, so `head(N)` is a slab of the embryo

In every released spatial stage, the row index follows one spatial axis almost
monotonically. Measured with `snippets/check_row_order.py` (5,000 cells; coverage is
the share of occupied 16³ voxels relative to a random 5,000-cell sample):

| File | Cells | max \|Spearman ρ\| of row index vs a principal axis | First 5,000 rows cover | Median contiguous 5,000-row block covers |
|---|---|---|---|---|
| E6.75 | 7,093 | 1.00 | 77 % | 77 % |
| E7.25 | 13,295 | 1.00 | 45 % | 47 % |
| E8.0 | 31,671 | 0.98 | 23 % | 25 % |
| E8.25_late | 58,716 | 0.99 | 10 % | 17 % |
| E8.75 | 24,826 | 0.92 | 21 % | 31 % |
| E9.5 | 53,742 | 0.82 | 9 % | 17 % |
| E9.5_mab21l2_ko | 50,294 | 0.92 | 15 % | 17 % |

The scRNA-seq files (`E8.5_RNA`, `E9.5_RNA`) have no coordinates, so this check does not apply to them.

So `adata[:n]`, `adata[start:start+n]`, or a loader with a row cap
(`read(..., max_rows=8000)`) does not give you "some cells of the stage". It
gives you one region of it. The file passes every format check, and the
spatial metrics then score it as a partial embryo. Two cases from our own
submissions, same stage, same method, only the sampling changed:

| Board | Cells submitted | Coverage vs a random sample | Real leaderboard score |
|---|---|---|---|
| T2 heart extrapolation | contiguous row block of E9.5 | 55 % | 38.6 |
| | random sample of E9.5 | 100 % | 48.94 |
| T2 embryo interpolation | first 8,000 rows of each stage (a row cap in the loader) | 61 % | 55.0 |
| | all rows, random sample | 100 % | 57.61 |

The row cap also biased the cell-type means used to build the prediction,
because they came from 8,000 of the later stage's 31,671 cells.

It is not a free win everywhere: on our T3 file a random sample scored lower
than the block, because one metric's validity gate flipped. Treat sampling as
something to check per board, and make it a deliberate choice instead of
an accident of `head()`.

```bash
python snippets/check_row_order.py data/raw/*.h5ad --n 5000
```

To see the region in 3D and get the same coverage number for your own
prediction file, use the
[Submission Viewer](https://github.com/Shashwat-srivastav/Virtual-Embryo-Submission-Viewer)
(`--compare my_prediction.h5ad source_stage.h5ad`).

## Bonus: unaudited precomputed embeddings

Some files ship precomputed embeddings (UMAP/PCA) alongside raw `X`. Their
exact provenance relative to any external-data rules in the Challenge isn't
documented, so treat them as untrusted and prefer recomputing from raw `.X`
yourself rather than consuming them directly.

## Once you're past these: a starter baseline

If you want a working, submission-ready starting point after handling the
gotchas above, see the companion contribution:
[Cell-Type-Conditioned Interpolation Baseline](https://github.com/Shashwat-srivastav/VirtualEmbryo-basline),
a cell-type-conditioned interpolation baseline, stronger than the naive
`copy_nearest`/`linear` starting points, with the same submission-schema
validation these snippets assume.

## Gotcha Doctor: paste an error, get the gotcha

`gotcha_doctor.py` (standard library only, no install, no network) matches an
error message or a symptom description against keyword lists for the seven
gotchas above and prints the most likely one with its fix:

```bash
python gotcha_doctor.py "adata_view.copy() raises ValueError: To copy an AnnData object in backed mode..."
python gotcha_doctor.py "I submitted the first 20000 rows and the shape metrics came out very low"
echo "obs_names has duplicates after concat" | python gotcha_doctor.py
```

```
Best match: backed_anndata_copy  (5 keyword hit(s): \.copy\(\), backed, to_memory, backed mode, pass a filename)
  an error when calling .copy() on an AnnData object opened with backed='r', or
  .to_memory() being unexpectedly slow

Suggested fix:
  .copy() isn't supported on a backed view; .to_memory() works but eagerly
  pulls all of obs/uns/varm too. Materialize only the fields you actually
  need (X, obsm) instead.
```

It is plain keyword matching (the `GOTCHAS` dict, about 60 lines of logic you
can read end to end), not a model. It shows which keywords fired so you can
judge the match yourself. A description that shares no vocabulary with any
gotcha gets "No match" instead of a guess, and a tie is reported as
a tie. It only knows these seven gotchas.

## Tests

```bash
python tests/test_gotcha_doctor.py
```

One realistic symptom per gotcha must map to the right one, a paraphrase with
no shared vocabulary must score zero, and the test cases must cover the whole
taxonomy.

## License

MIT. See [LICENSE](LICENSE).

---

Filed as a Community Contribution to the NeurIPS 2026 Virtual Embryo Challenge.
