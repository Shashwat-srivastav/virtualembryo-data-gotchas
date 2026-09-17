"""Verify an .h5ad file's var_names match the official gene panel byte-for-byte,
in order. Silent mismatches here mean every downstream gene-indexed array is
misaligned without erroring.

Usage:
    python check_gene_panel.py path/to/file.h5ad path/to/panel.genes.txt
"""
import sys

import anndata as ad


def check(h5ad_path: str, genes_file: str) -> bool:
    panel = [l.strip() for l in open(genes_file) if l.strip()]
    a = ad.read_h5ad(h5ad_path, backed="r")
    file_genes = list(a.var_names)

    if file_genes == panel:
        print(f"OK: {h5ad_path} var_names match {genes_file} exactly ({len(panel)} genes, same order)")
        return True

    print(f"MISMATCH: {h5ad_path} ({len(file_genes)} genes) vs {genes_file} ({len(panel)} genes)")
    missing = [g for g in panel if g not in file_genes]
    extra = [g for g in file_genes if g not in panel]
    if missing:
        print(f"  panel genes NOT in file (file wasn't measured on these): {missing[:10]}"
              f"{' ...' if len(missing) > 10 else ''}")
    if extra:
        print(f"  file genes NOT in panel: {extra[:10]}{' ...' if len(extra) > 10 else ''}")
    if not missing and not extra and file_genes != panel:
        print("  same gene SET but different ORDER -- reindex before using positionally")
    return False


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    ok = check(sys.argv[1], sys.argv[2])
    raise SystemExit(0 if ok else 1)
