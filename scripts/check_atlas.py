#!/usr/bin/env python3
"""Refuse an atlas that is not fit to publish.

`verify_genes.py --migrate` on its own rebuilds data/atlas_v3.json as an
unverified skeleton: every entry goes back to 'unverified'. That is the right
first step of the pipeline and the wrong thing to commit. This check catches
it, and also a run cut short by a network failure, which looks the same: a
sudden collapse in entries or in verified count.

Used by both workflows and by scripts/rebuild.py. Exit code 1 means do not
publish this file.
"""
import json
import sys

MIN_ENTRIES = 30
MIN_VERIFIED = 10


def problems(path: str = "data/atlas_v3.json") -> list[str]:
    a = json.load(open(path))
    entries = a.get("entries", [])
    s, n = a.get("stats", {}), len(entries)
    print(f"entries={n} verified={s.get('verified')} "
          f"partial={s.get('partial')} unverified={s.get('unverified')}")
    out = []
    if n < MIN_ENTRIES:
        out.append(f"only {n} entries; expected at least {MIN_ENTRIES}")
    if (s.get("verified") or 0) < MIN_VERIFIED:
        out.append(f"only {s.get('verified')} verified; expected at least {MIN_VERIFIED}")
    # Count from the entries too: stats can be stale if a step was skipped.
    counted = sum(1 for e in entries if (e.get("verification") or {}).get("status") == "verified")
    if counted < MIN_VERIFIED:
        out.append(f"only {counted} entries carry status 'verified' — this looks like a "
                   "--migrate-only skeleton. Run: python3 scripts/rebuild.py")
    for e in entries:
        if not e.get("id") or e.get("genes") is None:
            out.append(f"malformed entry: {e.get('id')!r}")
    return out


if __name__ == "__main__":
    p = problems(*sys.argv[1:])
    if p:
        print("NOT FIT TO PUBLISH:", *p, sep="\n  - ")
        sys.exit(1)
    print("ok")
