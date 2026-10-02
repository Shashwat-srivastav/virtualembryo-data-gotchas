"""Regression test for the keyword-matching classifier: one case per known
gotcha (using realistic phrasing that shares vocabulary with its category),
plus a paraphrase stress test with no shared keywords to confirm it
reports no match instead of guessing.
Run with: python test_gotcha_doctor.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gotcha_doctor import diagnose, GOTCHAS  # noqa: E402

CASES = {
    "gene_panel_mismatch": "My var_names has 500 genes but the official genes.txt panel "
                            "for this board only lists 498 -- two extra columns I didn't expect.",
    "duplicate_obs_names": "I grouped by obs_names and ended up with fewer unique cells "
                            "than rows in the file -- turns out some barcodes repeat.",
    "coordinate_frame_drift": "When I plot spatial_3D from two different stage files on "
                               "the same axes, one is centered near zero and the other is "
                               "way off at x=-5000, totally different bounding boxes.",
    "backed_anndata_copy": "adata_view.copy() raises ValueError: To copy an AnnData object "
                            "in backed mode, pass a filename, use .to_memory() instead.",
    "csc_slicing_cost": "Randomly subsampling 2000 rows from this h5ad file with backed='r' "
                         "takes forever, the sparse slicing is so slow.",
    "unharmonized_celltypes": "Endocardium shows up at E9.5 but there's no 'Endocardium' "
                               "label at all in the later stage's obs['celltype'] column.",
    "spatially_ordered_rows": "I submitted the first 20000 rows of E9.5 and the shape metric "
                               "came out very low, the cells look like one region in 3D.",
}

n_correct = 0
for expected, symptom in CASES.items():
    name, score, runner_up, hits = diagnose(symptom)
    ok = name == expected and score > 0
    n_correct += ok
    print(f"{'OK ' if ok else 'MISS'} expected={expected:24s} got={name:24s} score={score} hits={hits}")
assert n_correct == len(CASES), f"only {n_correct}/{len(CASES)} matched correctly"
print(f"{n_correct}/{len(CASES)} correct\n")

# A paraphrase sharing no vocabulary with any category should score 0
# everywhere instead of picking a wrong category.
name, score, runner_up, hits = diagnose(
    "my model behaves oddly when I merge two stage files together for training")
assert score == 0, f"expected zero keyword hits on a true paraphrase, got score={score} ({name})"
print(f"paraphrase-with-no-shared-vocabulary correctly scores 0 (would report 'No match'): PASS")

assert set(GOTCHAS) == set(CASES), "GOTCHAS taxonomy and test cases have drifted apart"
print("taxonomy coverage check: PASS")

print("\nALL TESTS PASS")
