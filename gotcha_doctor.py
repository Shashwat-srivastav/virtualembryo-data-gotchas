"""Paste an error message or a description of weird behavior from the Virtual
Embryo Challenge's .h5ad files, and get back the most likely known data
gotcha it matches, plus a pointer to the fix. The gotchas are the seven
documented in this repo's README.

Keyword matching with the standard library only: no install, no network
call. It will miss a description that avoids every listed keyword. See the
"Gotcha Doctor" section of the README for what it does and doesn't catch.

Usage:
    python gotcha_doctor.py "your error message or symptom description here"
    echo "your error message" | python gotcha_doctor.py
"""
import re
import sys

GOTCHAS = {
    "gene_panel_mismatch": {
        "criterion": "var_names don't match the official gene panel for the board "
                     "(wrong count, wrong order, or a gene the file never measured)",
        "fix": "Verify var_names against the official *.genes.txt panel byte-for-byte, "
               "in order, before treating columns as aligned to anything else.",
        "keywords": [
            r"var_names", r"gene panel", r"genes\.txt", r"gene (count|order)",
            r"\bpanel\b", r"gene.{0,15}(mismatch|missing|extra)", r"wrong (number|count) of genes",
            r"column.{0,10}(align|order|mismatch)",
        ],
    },
    "duplicate_obs_names": {
        "criterion": "obs_names / cell barcodes are not unique within a file, or a "
                     "join/index on obs_names merged unrelated cells",
        "fix": "Check for duplicate obs_names before joining or indexing on them; use "
               "positional indices or another unique key instead.",
        "keywords": [
            r"obs_names", r"\bbarcodes?\b", r"duplicate", r"not unique", r"non-unique",
            r"sample_id", r"repeat(s|ed)? (barcode|cell|row)",
        ],
    },
    "coordinate_frame_drift": {
        "criterion": "spatial_3D coordinates look wrong when pooled across files, or a "
                     "model trained on pooled coordinates behaves oddly",
        "fix": "Coordinates are per-file local, differing in both scale and origin. "
               "Re-center (and optionally rescale) each file to its own centroid "
               "before pooling coordinates across files.",
        "keywords": [
            r"spatial_3d", r"coordinate", r"centroid", r"bounding box", r"\borigin\b",
            r"\bscale\b", r"pool(ed|ing)? (coordinates|across)", r"registration",
            r"different (scale|offset|box)",
        ],
    },
    "backed_anndata_copy": {
        "criterion": "an error when calling .copy() on an AnnData object opened with "
                     "backed='r', or .to_memory() being unexpectedly slow",
        "fix": ".copy() isn't supported on a backed view; .to_memory() works but "
               "eagerly pulls all of obs/uns/varm too. Materialize only the fields "
               "you actually need (X, obsm) instead.",
        "keywords": [
            r"\.copy\(\)", r"backed", r"to_memory", r"backed mode", r"AnnData object",
            r"pass a filename",
        ],
    },
    "csc_slicing_cost": {
        "criterion": "slicing or subsetting rows from a large .h5ad file is much "
                     "slower than expected",
        "fix": "X is stored CSC (compressed-sparse-column), so scattered row indices are "
               "slow. Sort the random indices before slicing (a[np.sort(idx)]), or read X "
               "into memory once and sample from that. Do not switch to a contiguous "
               "block for anything you submit or score: rows are in spatial order "
               "(see spatially_ordered_rows).",
        "keywords": [
            r"\bcsc\b", r"sparse", r"slic(e|ing)", r"subsampl",
            r"slow(er|ly)?", r"takes forever", r"random.{0,10}(index|subsample|rows?)",
        ],
    },
    "unharmonized_celltypes": {
        "criterion": "a cell-type label present at one stage appears 'missing' at "
                     "another stage, or a groupby/diff on celltype gives surprising results",
        "fix": "Cell-type vocabulary isn't harmonized across stages, so a missing label "
               "may be a rename, not a biological absence. Diff the label sets between "
               "stages and eyeball for near-matches before trusting the absence.",
        "keywords": [
            r"cell.?type", r"celltype", r"groupby", r"label.{0,15}(missing|absent)",
            r"vocabulary", r"rename", r"no .* label",
        ],
    },
    "spatially_ordered_rows": {
        "criterion": "a prediction built from the first N rows, a contiguous row block, or a "
                     "row-capped loader scores low on the spatial or shape metrics, or its "
                     "cells look like one region of the embryo",
        "fix": "Rows are stored in spatial order, so head(N) or a[start:start+n] is a slab "
               "of the embryo. Take a random sample of all rows instead, remove row caps "
               "from loaders, and check coverage with snippets/check_row_order.py or the "
               "Submission Viewer's --compare.",
        "keywords": [
            r"first \d[\d,]* (rows|cells)", r"head\(", r"\[:\s*\d", r"row block", r"contiguous",
            r"row cap", r"max_rows", r"slab", r"one region", r"partial embryo", r"coverage",
            r"shape metric", r"d2_shape|occupancy|\bsdd\b|\bods\b", r"spatial (metric|score)s?.{0,20}low",
        ],
    },
}


def diagnose(symptom: str):
    """Returns (best_name, score, runner_up_score, matched_keywords)."""
    text = symptom.lower()
    scores = {}
    matches = {}
    for name, info in GOTCHAS.items():
        hits = [kw for kw in info["keywords"] if re.search(kw, text)]
        scores[name] = len(hits)
        matches[name] = hits

    ranked = sorted(scores.items(), key=lambda kv: -kv[1])
    best_name, best_score = ranked[0]
    runner_up_score = ranked[1][1] if len(ranked) > 1 else 0
    return best_name, best_score, runner_up_score, matches[best_name]


def main():
    if len(sys.argv) > 1:
        symptom = " ".join(sys.argv[1:])
    else:
        symptom = sys.stdin.read().strip()
    if not symptom:
        raise SystemExit(__doc__)

    best_name, best_score, runner_up, hit_keywords = diagnose(symptom)

    if best_score == 0:
        print("No match: none of the seven known gotchas' keywords appear in this text.")
        print("See the companion virtualembryo-data-gotchas repo for the full list, "
              "or this may be a different issue entirely.")
        return

    match = GOTCHAS[best_name]
    print(f"Best match: {best_name}  ({best_score} keyword hit(s): {', '.join(hit_keywords)})")
    print(f"  {match['criterion']}")
    print(f"\nSuggested fix:\n  {match['fix']}")
    if best_score == runner_up:
        print(f"\nNote: tied with another category on keyword count, so this is a rough "
              f"match. Read the description above and judge for yourself.")


if __name__ == "__main__":
    main()
