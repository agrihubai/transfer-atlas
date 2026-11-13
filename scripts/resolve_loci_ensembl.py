#!/usr/bin/env python3
"""
Resolve genome-project locus tags against Ensembl Plants annotation files.

Why this exists
---------------
Locus tags like Csa5G623470, MELO3C014437 or Cla97C02G036390 are not indexed in
NCBI and Sol Genomics only covers Solanaceae. Ensembl Plants has the rest, but
its REST service is gone: rest.ensemblgenomes.org serves a certificate for
another host, and plants.ensembl.org now redirects to www.ensembl.org. What does
work is the annotation itself, over HTTPS from ftp.ebi.ac.uk.

So this script streams the GFF3 for each species, pulls out only the gene lines
for the tags asked about, and writes the answers to data/locus_resolutions.json.
Nothing large is kept: the resolutions file is small, committed, and carries the
build and the URL it came from, so the result is auditable without re-downloading
a genome and the verification pipeline needs no network for it.

Usage
-----
  python3 scripts/resolve_loci_ensembl.py            # resolve what the atlas needs
  python3 scripts/resolve_loci_ensembl.py --list     # show the pending tags only
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = DATA / "locus_resolutions.json"
FTP = "https://ftp.ebi.ac.uk/ensemblgenomes/pub/current/plants/gff3"
UA = "TransferAtlas/3.0 (https://agrihubai.github.io/transfer-atlas)"

# species directory -> annotation file. Ensembl's own names; update when the
# release changes (the version number is part of the filename).
SPECIES = {
    "Cucumis sativus": ("cucumis_sativus", "Cucumis_sativus.ASM407v2.63.chr.gff3.gz"),
    "Cucumis melo": ("cucumis_melo", "Cucumis_melo.Melonv4.63.gff3.gz"),
    "Citrullus lanatus": ("citrullus_lanatus", "Citrullus_lanatus.Cla97_v1.63.chr.gff3.gz"),
    "Brassica rapa": ("brassica_rapa", "Brassica_rapa.Brapa_1.0.63.chr.gff3.gz"),
    "Brassica oleracea": ("brassica_oleracea", None),
    "Capsicum annuum": ("capsicum_annuum", None),
    "Daucus carota": ("daucus_carota", None),
}


def variants(locus: str) -> list[str]:
    """Spellings of the same tag across builds.

    Ensembl writes cucumber genes as Csa_5G623470 where the literature and the
    atlas write Csa5G623470. Without this the tags simply miss.
    """
    out = [locus]
    m = re.match(r"^(Csa)(\d+G\d+)$", locus, re.I)
    if m:
        out.append(f"{m.group(1)}_{m.group(2)}")
    out.append(re.sub(r"\.\d+$", "", locus))
    seen, uniq = set(), []
    for v in out:
        if v.lower() not in seen:
            seen.add(v.lower())
            uniq.append(v)
    return uniq


def gene_lines(species_dir: str, filename: str, wanted: dict[str, str]):
    """Stream one annotation file, yielding gene rows that match a wanted tag."""
    url = f"{FTP}/{species_dir}/{filename}"
    pattern = re.compile("|".join(re.escape(v) for v in wanted), re.I)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=300) as r:
        with gzip.GzipFile(fileobj=io.BufferedReader(r)) as gz:  # type: ignore[arg-type]
            for raw in gz:
                line = raw.decode("utf-8", "replace")
                if "\tgene\t" not in line or not pattern.search(line):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                attrs = dict(
                    (kv.split("=", 1) + [""])[:2] for kv in f[8].split(";") if "=" in kv
                )
                yield {
                    "ensembl_id": (attrs.get("gene_id") or attrs.get("ID", "")).replace("gene:", ""),
                    "name": attrs.get("Name", ""),
                    "description": attrs.get("description", ""),
                    "biotype": attrs.get("biotype", ""),
                    "chromosome": f[0],
                    "start": f[3],
                    "end": f[4],
                    "source_url": url,
                }


def pending() -> dict[str, list[str]]:
    """Locus tags in the atlas that nothing has resolved yet, grouped by species."""
    atlas = json.loads((DATA / "atlas_v3.json").read_text(encoding="utf-8"))
    out: dict[str, list[str]] = {}
    for e in atlas.get("entries", []):
        sp = e.get("scientific_name") or ""
        for r in (e.get("locus_check") or {}).get("results", []):
            if r.get("state") == "confirmed_in_species":
                continue
            out.setdefault(sp, [])
            if r["locus"] not in out[sp]:
                out[sp].append(r["locus"])
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--list", action="store_true", help="show pending tags and exit")
    a = ap.parse_args()

    todo = pending()
    if a.list:
        for sp, loci in sorted(todo.items()):
            print(f"{sp:<24}{', '.join(loci)}")
        return 0

    resolved = {}
    if OUT.exists():
        resolved = json.loads(OUT.read_text(encoding="utf-8")).get("loci", {})
    today = time.strftime("%Y-%m-%d")

    for sp, loci in sorted(todo.items()):
        entry = SPECIES.get(sp)
        if not entry or not entry[1]:
            print(f"-- {sp}: no annotation file configured; skipped ({', '.join(loci)})")
            continue
        sdir, fname = entry
        wanted = {v: locus for locus in loci for v in variants(locus)}
        print(f"-- {sp}: streaming {fname} for {len(loci)} tag(s)")
        found = 0
        try:
            for row in gene_lines(sdir, fname, wanted):
                key = next((wanted[v] for v in wanted
                            if v.lower() in row["ensembl_id"].lower()), None)
                if not key:
                    continue
                row.update({"species": sp, "build": fname.split(".")[1],
                            "checked_on": today, "source": "Ensembl Plants annotation"})
                resolved[key] = row
                found += 1
                print(f"   {key:<20}{row['ensembl_id']:<20}chr{row['chromosome']:<6}"
                      f"{row['name'] or '(no gene name)'}")
        except Exception as e:  # network or gzip trouble: report, do not crash
            print(f"   ! {type(e).__name__}: {e}", file=sys.stderr)
        if found < len(loci):
            missing = [l for l in loci if l not in resolved]
            if missing:
                print(f"   not in this build: {', '.join(missing)}")

    OUT.write_text(json.dumps({
        "_about": ("Locus tags resolved against Ensembl Plants annotation by "
                   "scripts/resolve_loci_ensembl.py. Each record carries the build and the "
                   "file it came from. A tag absent here was not in that build — which is "
                   "not proof the tag is wrong, only that this build does not carry it."),
        "loci": resolved,
    }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"\nwrote {OUT.relative_to(ROOT)} with {len(resolved)} resolution(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
