#!/usr/bin/env python3
"""Patent/IP signal layer for the Transfer Atlas.

For every gene in data/atlas.json this always emits a live Google Patents
search link (no key needed) so a breeder can see the IP landscape in one click.
If a PatentsView API key is present (secret PATENTSVIEW_API_KEY), it also
pulls an approximate US patent count and sample titles as an automated signal.

Honesty: a count or a link is a research signal, never a freedom-to-operate
opinion. "Patent activity exists" is not "this use is barred".
"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results"
KEY = os.environ.get("PATENTSVIEW_API_KEY", "")
PV = "https://search.patentsview.org/api/v1/patent/"


def gene_term(gene):
    return gene.split(" (")[0].strip()


def trait_word(trait):
    return {"salt": "salt tolerance", "seedless": "seedless"}.get(trait, trait)


def google_link(gene, trait):
    q = gene_term(gene) + " " + trait_word(trait)
    return "https://patents.google.com/?q=" + urllib.parse.quote(q)


def patentsview_count(term):
    """Best-effort US patent count via PatentsView. Returns (count, sample) or (None, '')."""
    if not KEY:
        return None, ""
    q = json.dumps({"_text_phrase": {"patent_abstract": term}})
    f = json.dumps(["patent_id", "patent_title"])
    o = json.dumps({"size": 5})
    url = PV + "?q=" + urllib.parse.quote(q) + "&f=" + urllib.parse.quote(f) + "&o=" + urllib.parse.quote(o)
    req = urllib.request.Request(url, headers={"X-Api-Key": KEY, "Accept": "application/json"})
    try:
        time.sleep(1.0)
        with urllib.request.urlopen(req, timeout=60) as r:
            d = json.loads(r.read())
        count = d.get("total_hits") or d.get("count") or (len(d.get("patents", [])) if d.get("patents") else 0)
        titles = [p.get("patent_title", "")[:70] for p in (d.get("patents") or [])[:3]]
        return count, " | ".join(titles)
    except Exception as e:  # noqa: BLE001
        print(f"  PatentsView error for '{term}': {e}", file=sys.stderr)
        return None, ""


def main():
    atlas = json.loads((ROOT / "data" / "atlas.json").read_text())
    OUT.mkdir(exist_ok=True)
    rows = []
    for e in atlas["entries"]:
        term = gene_term(e["gene"])
        link = google_link(e["gene"], e["trait"])
        count, sample = patentsview_count(term)
        rows.append({
            "gene": e["gene"], "trait": e["trait"],
            "patent_count": "n/a" if count is None else count,
            "source": "manual-link" if count is None else "patentsview",
            "search_link": link, "sample_titles": sample,
        })
        print(f"  {term:<16} count={rows[-1]['patent_count']}  {link}")

    out = OUT / "patents.tsv"
    import csv
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene", "trait", "patent_count", "source", "search_link", "sample_titles"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    mode = "PatentsView API + links" if KEY else "Google Patents links only (no API key set)"
    print(f"\nWrote {len(rows)} rows -> {out.name}  [{mode}]")


if __name__ == "__main__":
    main()
