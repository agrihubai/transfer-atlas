#!/usr/bin/env python3
"""
Resolve genome-project locus tags against public genome annotation files.

Why this exists
---------------
Locus tags like Csa5G623470, MELO3C014437 or Cla97C02G036390 are not indexed in
NCBI and Sol Genomics only covers Solanaceae. Ensembl Plants has the rest, but
its REST service is gone: rest.ensemblgenomes.org serves a certificate for
another host, and plants.ensembl.org now redirects to www.ensembl.org. What does
work is the annotation itself, over HTTPS: Ensembl Genomes' files on ftp.ebi.ac.uk
for cucurbits and Brassicas, and Sol Genomics' own FTP tree for eggplant, whose
tags are in no other reachable database.

So this script streams the GFF3 for each species, pulls out only the gene lines
for the tags asked about, and writes the answers to data/locus_resolutions.json.
Nothing large is kept: the resolutions file is small, committed, and carries the
build and the URL it came from, so the result is auditable without re-downloading
a genome and the verification pipeline needs no network for it.

Usage
-----
  python3 scripts/resolve_loci.py            # resolve what the atlas needs
  python3 scripts/resolve_loci.py --list     # show the pending tags only
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = DATA / "locus_resolutions.json"
FTP = "https://ftp.ebi.ac.uk/ensemblgenomes/pub/current/plants/gff3"
UA = "TransferAtlas/3.0 (https://agrihubai.github.io/transfer-atlas)"

# Where each species' annotation lives. Ensembl's filenames carry the release
# number, so they change when the release does. gzip is inferred from the suffix.
FTP = "https://ftp.ebi.ac.uk/ensemblgenomes/pub/current/plants/gff3"
SGN_FTP = "https://solgenomics.net/ftp/genomes"

SPECIES = {
    "Cucumis sativus": f"{FTP}/cucumis_sativus/Cucumis_sativus.ASM407v2.63.chr.gff3.gz",
    "Cucumis melo": f"{FTP}/cucumis_melo/Cucumis_melo.Melonv4.63.gff3.gz",
    "Citrullus lanatus": f"{FTP}/citrullus_lanatus/Citrullus_lanatus.Cla97_v1.63.chr.gff3.gz",
    "Brassica rapa": f"{FTP}/brassica_rapa/Brassica_rapa.Brapa_1.0.63.chr.gff3.gz",
    # Eggplant is absent from Ensembl Plants; Sol Genomics hosts the V4 annotation,
    # and this file carries functional descriptions as well as gene models.
    "Solanum melongena": (f"{SGN_FTP}/Solanum_melongena_V4_Pangenome/Annotation_V4/"
                          "SMEL_V4_function_IPR.gff3"),
}


def variants(locus: str) -> list[str]:
    """Spellings of the same tag across builds and sources.

    Three differences seen in practice: Ensembl writes cucumber genes as
    Csa_5G623470 where the literature writes Csa5G623470; the eggplant V4 files
    write SMEL4_08g009050 where the atlas writes SMEL4.1_08g009050; and Sol
    Genomics appends a model suffix (.1). Without these the tags simply miss.
    """
    out = [locus]
    m = re.match(r"^(Csa)(\d+G\d+)$", locus, re.I)
    if m:
        out.append(f"{m.group(1)}_{m.group(2)}")
    if locus.upper().startswith("SMEL4.1_"):
        out.append("SMEL4_" + locus.split("_", 1)[1])
    base = re.sub(r"\.\d+$", "", locus)
    out += [base, base + ".1"]
    seen, uniq = set(), []
    for v in out:
        if v and v.lower() not in seen:
            seen.add(v.lower())
            uniq.append(v)
    return uniq


def gene_lines(url: str, wanted: dict[str, str]):
    """Stream one annotation file, yielding gene rows that match a wanted tag."""
    pattern = re.compile("|".join(re.escape(v) for v in wanted), re.I)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=300) as r:
        stream = (gzip.GzipFile(fileobj=io.BufferedReader(r))  # type: ignore[arg-type]
                  if url.endswith(".gz") else r)
        for raw in stream:
            line = raw.decode("utf-8", "replace")
            if "\tgene\t" not in line or not pattern.search(line):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attrs = dict(
                (kv.split("=", 1) + [""])[:2] for kv in f[8].split(";") if "=" in kv
            )
            gid = (attrs.get("gene_id") or attrs.get("ID", "")).replace("gene:", "")
            # Ensembl puts the function in description=, Sol Genomics in Note=
            desc = urllib.parse.unquote(attrs.get("description") or attrs.get("Note") or "")
            name = urllib.parse.unquote(attrs.get("Name", ""))
            # Sol Genomics sets Name= to the gene id itself; a name that only
            # repeats the identifier carries no information, so drop it.
            if name and name.split(".")[0].lower() == gid.split(".")[0].lower():
                name = ""
            yield {
                "ensembl_id": urllib.parse.unquote(gid),
                "name": name,
                "description": desc,
                "biotype": attrs.get("biotype", ""),
                "chromosome": f[0],
                "start": f[3],
                "end": f[4],
                "source_url": url,
            }


def pending() -> dict[str, list[str]]:
    """Locus tags the atlas needs and nothing has resolved yet, grouped by species.

    Reads three places, because any one of them can be empty: the locus checks in
    atlas_v3.json after a verification run, the raw gene_id claims (so this works
    on a freshly migrated atlas too), and the gene_locus maps in corrections.json,
    where a curator has pinned a symbol to a tag the draft never carried.
    """
    atlas = json.loads((DATA / "atlas_v3.json").read_text(encoding="utf-8"))
    try:
        fixes = json.loads((DATA / "corrections.json").read_text(encoding="utf-8"))["entries"]
    except (OSError, KeyError, json.JSONDecodeError):
        fixes = {}
    done = set()
    if OUT.exists():
        try:
            done = set(json.loads(OUT.read_text(encoding="utf-8")).get("loci", {}))
        except json.JSONDecodeError:
            pass

    out: dict[str, list[str]] = {}

    def add(sp: str, locus: str) -> None:
        if not sp or not locus or locus in done:
            return
        out.setdefault(sp, [])
        if locus not in out[sp]:
            out[sp].append(locus)

    for e in atlas.get("entries", []):
        sp = e.get("scientific_name") or ""
        checks = (e.get("locus_check") or {}).get("results", [])
        if checks:
            for r in checks:
                if r.get("state") != "confirmed_in_species":
                    add(sp, r["locus"])
        else:
            for locus in re.findall(r"[A-Za-z]{3,}[0-9][A-Za-z0-9._\-]{4,}",
                                    (e.get("claims") or {}).get("gene_id", "") or ""):
                add(sp, locus)
        for locus in ((fixes.get(e.get("id")) or {}).get("gene_locus") or {}).values():
            add(sp, locus)
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
        url = SPECIES.get(sp)
        if not url:
            print(f"-- {sp}: no annotation file configured; skipped ({', '.join(loci)})")
            continue
        wanted = {v: locus for locus in loci for v in variants(locus)}
        print(f"-- {sp}: streaming {url.rsplit('/', 1)[-1]} for {len(loci)} tag(s)")
        found = 0
        try:
            for row in gene_lines(url, wanted):
                key = next((wanted[v] for v in wanted
                            if v.lower() in row["ensembl_id"].lower()), None)
                if not key:
                    continue
                row.update({"species": sp, "build": url.rsplit("/", 1)[-1],
                            "checked_on": today, "source": "Ensembl Plants annotation"})
                resolved[key] = row
                found += 1
                print(f"   {key:<20}{row['ensembl_id']:<20}chr{row['chromosome']:<6}"
                      f"{row['name'] or row['description'][:52] or '(no annotation)'}")
        except Exception as e:  # network or gzip trouble: report, do not crash
            print(f"   ! {type(e).__name__}: {e}", file=sys.stderr)
        if found < len(loci):
            missing = [l for l in loci if l not in resolved]
            if missing:
                print(f"   not in this build: {', '.join(missing)}")

    OUT.write_text(json.dumps({
        "_about": ("Locus tags resolved against public genome annotation files by "
                   "scripts/resolve_loci.py. Each record carries the build and the "
                   "file it came from. A tag absent here was not in that build — which is "
                   "not proof the tag is wrong, only that this build does not carry it."),
        "loci": resolved,
    }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"\nwrote {OUT.relative_to(ROOT)} with {len(resolved)} resolution(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
