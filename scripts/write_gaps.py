#!/usr/bin/env python3
"""
Write docs/verification-gaps.md from the current atlas state.

A hand-written gap list goes stale the moment the pipeline runs again, and a
stale list is worse than none: it sends people looking for things that are
already resolved. So the list is generated, ordered by what a second opinion can
actually close.

  python3 scripts/write_gaps.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ATLAS = ROOT / "data/atlas_v3.json"
OUT = ROOT / "docs/verification-gaps.md"

CONFIRMED_MATCH = ("author+year+gene", "author+year+journal")


def gene_unresolved(g: dict) -> bool:
    v = g.get("verification") or {}
    return not ((v.get("ncbi_gene") or {}).get("found") or v.get("genbank_locus"))


def confirmed_refs(e: dict) -> list[dict]:
    out = []
    for g in e.get("genes", []):
        pv = (g.get("verification") or {}).get("pubmed") or {}
        if pv.get("found") and str(pv.get("match_type", "")).startswith(CONFIRMED_MATCH):
            out.append(pv["hits"][0])
    for c in (e.get("verification") or {}).get("curated_references") or []:
        if c.get("pubmed_title"):
            out.append({"pmid": c["pmid"], "title": c["pubmed_title"],
                        "journal": c.get("pubmed_journal", "")})
    return out


def awaiting(e: dict) -> list[dict]:
    out = []
    for g in e.get("genes", []):
        pv = (g.get("verification") or {}).get("pubmed") or {}
        if str(pv.get("match_type", "")).startswith("awaiting"):
            out.append(pv["hits"][0])
    return out


def main() -> int:
    atlas = json.loads(ATLAS.read_text(encoding="utf-8"))
    entries = atlas["entries"]
    s = atlas.get("stats", {})
    L: list[str] = []
    A = L.append

    A("# What is still open, and what would close it\n\n")
    A(f"Generated from `data/atlas_v3.json` on {time.strftime('%d %B %Y')} by "
      "`scripts/write_gaps.py`. Do not edit by hand — rerun the script.\n\n")
    A(f"State: **{s.get('verified', 0)} verified, {s.get('partial', 0)} partial, "
      f"{s.get('unverified', 0)} unverified** of {s.get('total_entries', 0)} entries.\n\n")

    A("## What counts, so nobody wastes effort\n\n")
    A("A reference counts only when PubMed matches the cited **author, year and journal**, or "
      "when a supplied PMID resolves. A PMID found by topic search alone does not count: it may "
      "be a different paper about the same gene.\n\n")
    A("**The most useful thing anyone can send** is a gene symbol, its crop species, and a DOI "
      "or PMID for a paper that functionally tests that gene in that crop. The pipeline resolves "
      "the rest. Reconstructed citations do not help — an author and year with no findable "
      "journal is what the original draft had, and it is what produced most of these gaps.\n\n")
    A("**Not genes, so do not look for them:** `uORF` (a 5' UTR feature), `ROS scavengers`, "
      "`ABA catabolism`, `prickle loci`. These were descriptive phrases in the draft's gene "
      "field that the parser read as symbols.\n\n")

    # 1. awaiting a human glance
    aw = [(e, awaiting(e)) for e in entries if awaiting(e)]
    if aw:
        A(f"## Needs one look, nothing more ({len(aw)})\n\n")
        A("The author and journal match the citation, but the title does not carry the trait, so "
          "the pipeline will not count it. A person can settle this in seconds.\n\n")
        for e, hits in aw:
            h = hits[0]
            A(f"- **{e['id']}** — {e['crop']}, {e['trait']}\n")
            A(f"  - cited: {e['claims']['scientific_reference']}\n")
            A(f"  - candidate: PMID {h['pmid']} ({h.get('journal','')}) {h['title'][:92]}\n")
            A("  - question: is this the intended paper? yes or no\n")

    # 2. confirmed paper, unresolved gene identity
    p1 = [e for e in entries
          if confirmed_refs(e) and all(gene_unresolved(g) for g in e["genes"])]
    if p1:
        A(f"\n## A confirmed paper exists, but no gene symbol resolves ({len(p1)})\n\n")
        A("The cheapest wins: the evidence is confirmed, only the identifier is missing. The "
          "paper's own data-availability section usually has it.\n\n")
        for e in p1:
            A(f"- **{e['id']}** — {e['crop']} ({e['scientific_name']}), {e['trait']}\n")
            A(f"  - unresolved: {', '.join('`' + g['symbol'] + '`' for g in e['genes'])}\n")
            r = confirmed_refs(e)[0]
            A(f"  - confirmed paper: PMID {r['pmid']} — {r['title'][:88]}\n")
            A("  - needed: an accession or gene ID for those symbols, read off a record\n")

    # 3. genes resolve, no confirmed reference
    p2 = [e for e in entries
          if not confirmed_refs(e) and any(not gene_unresolved(g) for g in e["genes"])]
    if p2:
        A(f"\n## Genes resolve, no reference is confirmed ({len(p2)})\n\n")
        A("Identity is settled; the functional evidence is not. A DOI or PMID for a knockout, "
          "overexpression, allele-swap or field study closes these.\n\n")
        for e in p2:
            A(f"- **{e['id']}** — {e['crop']}, {e['trait']} (`{e['gene_label']}`)\n")
            A(f"  - cited: \"{e['claims']['scientific_reference']}\" — not matched\n")

    # 4. nothing at all
    p3 = [e for e in entries if not confirmed_refs(e) and all(gene_unresolved(g) for g in e["genes"])]
    if p3:
        A(f"\n## Nothing resolves: neither gene nor reference ({len(p3)})\n\n")
        A("These should be filled with a real source or retired from the atlas.\n\n")
        for e in p3:
            A(f"- **{e['id']}** — {e['crop']}, {e['trait']} (`{e['gene_label']}`)\n")
            A(f"  - cited: \"{e['claims']['scientific_reference']}\"\n")

    # 5. partially resolved multi-gene entries
    p4 = [e for e in entries
          if any(gene_unresolved(g) for g in e["genes"])
          and any(not gene_unresolved(g) for g in e["genes"])]
    if p4:
        A(f"\n## Some genes of an entry resolve, others do not ({len(p4)})\n\n")
        for e in p4:
            miss = [g["symbol"] for g in e["genes"] if gene_unresolved(g)]
            A(f"- **{e['id']}** — {e['crop']} ({e['verification']['genes_resolved']} resolved): "
              f"{', '.join('`' + m + '`' for m in miss)}\n")

    # 6. locus tags
    bad = [(e, r) for e in entries
           for r in (e.get("locus_check") or {}).get("results", [])
           if r.get("state") != "confirmed_in_species"]
    if bad:
        A(f"\n## Locus tags that no reachable database carries ({len(bad)})\n\n")
        A("Not evidence that a tag is wrong — the database that holds it is unreachable, or the "
          "tag comes from a different genome build than the one available.\n\n")
        for e, r in bad:
            A(f"- `{r['locus']}` — {e['id']}, {e['scientific_name']}\n")

    # 7. commercial status still open
    comm = [e for e in entries
            if (e.get("correction") or {}).get("commercial_status")]
    if comm:
        A(f"\n## Commercial status, which only a person can settle ({len(comm)})\n\n")
        A("Level 1 requires a named cultivar or product whose trait comes from **this entry's "
          "genetic solution**, not the same trait by another route.\n\n")
        for e in comm:
            A(f"- **{e['id']}** — {e['crop']}, {e['trait']}\n")
            A(f"  - {(e['correction']['commercial_status'])[:240]}\n")

    OUT.write_text("".join(L), encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)} ({len(''.join(L))} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
