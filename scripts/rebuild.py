#!/usr/bin/env python3
"""Regenerate data/atlas_v3.json the only safe way: the whole pipeline.

Runs, in order, the same steps as .github/workflows/verify.yml:
  verify_genes.py --migrate  ->  resolve_loci.py  ->  verify_genes.py --verify
  ->  verify_genes.py --accessions  ->  write_gaps.py  ->  check_atlas.py

Use this after editing data/corrections.json or any other curated input.
Running --migrate alone leaves every entry 'unverified'; committing that file
empties the published atlas.

Database answers are cached, so a re-run after a small edit is quick.
"""
import subprocess
import sys

STEPS = [
    ["scripts/verify_genes.py", "--migrate"],
    ["scripts/resolve_loci.py"],
    ["scripts/verify_genes.py", "--verify"],
    ["scripts/verify_genes.py", "--accessions"],
    ["scripts/write_gaps.py"],
    ["scripts/check_atlas.py"],
]

for step in STEPS:
    print(f"\n== {' '.join(step)}", flush=True)
    r = subprocess.run([sys.executable, *step])
    if r.returncode != 0:
        print(f"\nstopped: {' '.join(step)} exited with {r.returncode}. "
              "data/atlas_v3.json is not fit to commit.", file=sys.stderr)
        sys.exit(r.returncode)
print("\nrebuild complete; data/atlas_v3.json is fit to commit.")
