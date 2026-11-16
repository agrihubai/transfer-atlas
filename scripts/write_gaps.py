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

    A("## Already searched — please do not repeat these\n\n")
    A("Every unresolved symbol below has been looked for, by script, in all of these:\n\n")
    A("- **NCBI Gene**, as written and with the species prefix stripped (`SlGAD3` and `GAD3`), "
      "constrained to the stated organism, with hits ranked by whether the symbol is the "
      "official name or only an alias.\n")
    A("- **NCBI nuccore and protein**, by locus tag, constrained by organism.\n")
    A("- **Sol Genomics Network** feature search, for Solanaceae symbols and locus tags.\n")
    A("- **Ensembl Plants annotation files** (ftp.ebi.ac.uk), streamed per species and searched "
      "both by locus tag and by gene name — cucumber, melon, watermelon, both Brassicas, pepper "
      "and carrot. Cucumber tags also tried in Ensembl's underscore form (`Csa_5G623470`).\n")
    A("- **Eggplant V4 annotation** (hosted by Sol Genomics), including the V3-to-V4 conversion "
      "table.\n\n")
    A("These annotations mostly carry no gene names, which is why a name search found almost "
      "nothing: of 21 unresolved symbols it resolved none. So the useful answer is **an "
      "accession or gene ID read off a record or a paper's data-availability section** — not "
      "another database name, and not a guess.\n\n")
    A("**What is unreachable from the build environment**, so do not suggest it: Ensembl's REST "
      "service for plants (rest.ensemblgenomes.org serves a certificate for another host), "
      "plants.ensembl.org (redirects to an unreachable host), cucurbitgenomics.org (no working "
      "HTTPS), and the watermelon 97103 v2 build (nowhere reachable).\n\n")
    A("**Three rules, from what went wrong last time.** Identifiers supplied earlier included a "
      "fungal ITS sequence given as a watermelon gene, a bacterial gene and a tomato protein "
      "given as eggplant polyphenol oxidases, an NCBI id for a protein kinase given as a lettuce "
      "GGP2, and two PMIDs that were a cellulose-microfibril paper and a cardiovascular review. "
      "The PMIDs supplied were otherwise reliable. So:\n\n")
    A("1. Never give an identifier not read off a record. Write **unknown** instead — an honest "
      "gap is useful, a wrong accession costs more than it saves.\n")
    A("2. Always name the **journal** with author and year. The matcher requires it, because "
      "author and year alone cannot separate two papers.\n")
    A("3. Say what a paper does **not** cover. That is what saved several entries from false "
      "claims.\n\n")

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

    # 5b. the unresolved symbols themselves, grouped by species
    import re as _re
    PRE = ("Smel", "Sl", "Ca", "Cs", "Cm", "Cl", "Sm", "Ls", "Bo", "Br", "Dc", "Ac", "Le")

    def _strip(sym: str) -> str:
        for pre in sorted(PRE, key=len, reverse=True):
            m = _re.match(r"^" + pre + r"([A-Z0-9].*)$", sym)
            if m and len(m.group(1)) >= 2:
                return m.group(1)
        return sym

    by_sp: dict[str, list[tuple[str, str, str]]] = {}
    for e in entries:
        for g in e["genes"]:
            if gene_unresolved(g):
                by_sp.setdefault(e.get("scientific_name") or "unknown species", []).append(
                    (e["id"], g["symbol"], _strip(g["symbol"])))
    if by_sp:
        total = sum(len(v) for v in by_sp.values())
        A(f"\n## Every gene symbol with no identifier, by species ({total})\n\n")
        A("For each: the entry, the symbol as the atlas writes it, and the name that was "
          "searched after stripping the species prefix. An accession, a gene ID, or `unknown`.\n")
        for sp in sorted(by_sp):
            A(f"\n**{sp}**\n\n")
            for eid, sym, base in by_sp[sp]:
                A(f"- `{sym}` ({eid}) — searched as `{base}`\n")

    # 6. locus tags
    bad = [(e, r) for e in entries
           for r in (e.get("locus_check") or {}).get("results", [])
           if r.get("state") != "confirmed_in_species"]
    if bad:
        A(f"\n## Locus tags that no reachable database carries ({len(bad)})\n\n")
        A("Not evidence that a tag is wrong — the database that holds it is unreachable, or the "
          "tag comes from a different genome build than the one available. **The useful answer "
          "here is which genome build each tag belongs to**, so it can be converted.\n\n")
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
