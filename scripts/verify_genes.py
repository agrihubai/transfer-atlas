#!/usr/bin/env python3
"""
Transfer Atlas — gene & reference verification pipeline.

What this script does
---------------------
For every gene in the atlas it asks public databases what they actually hold,
and records the answer. It NEVER writes a gene ID, PMID, DOI or function that
it did not retrieve from a database. Anything it cannot confirm is left empty
and the entry is marked `unverified`.

Sources queried
---------------
  NCBI Gene      (E-utilities esearch/esummary)  -> gene ID, official symbol, description
  NCBI PubMed    (E-utilities esearch/esummary)  -> PMID + title + DOI for the cited reference
  Ensembl Plants (REST lookup/id, lookup/symbol) -> stable gene ID, biotype, coordinates

Usage
-----
  python3 scripts/verify_genes.py --migrate          # build v3 skeleton from v2 (offline)
  python3 scripts/verify_genes.py --verify            # query databases, fill verification blocks
  python3 scripts/verify_genes.py --verify --only Tomato
  python3 scripts/verify_genes.py --report            # print the summary table only

Network
-------
Needs outbound HTTPS to eutils.ncbi.nlm.nih.gov and rest.ensembl.org.
Responses are cached in .cache/ so re-runs are cheap and offline-friendly.
Set NCBI_API_KEY in the environment to raise the rate limit from 3 to 10 req/s.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
CACHE = ROOT / ".cache"
V2 = DATA / "atlas_v2.json"
V3 = DATA / "atlas_v3.json"
REPORT = DATA / "verification_report.tsv"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENSEMBL = "https://rest.ensembl.org"
UA = "TransferAtlas/3.0 (https://agrihubai.github.io/transfer-atlas; verification pipeline)"

# Scientific name per crop label, used to constrain NCBI Gene searches.
SPECIES = {
    "Tomato": "Solanum lycopersicum",
    "Pepper": "Capsicum annuum",
    "Cucumber": "Cucumis sativus",
    "Melon": "Cucumis melo",
    "Watermelon": "Citrullus lanatus",
    "Eggplant": "Solanum melongena",
    "Lettuce": "Lactuca sativa",
    "Brassica (Broccoli/Cabbage/Cauliflower)": "Brassica oleracea",
    "Carrot": "Daucus carota",
    "Onion": "Allium cepa",
}

# Ensembl Plants species keys for symbol lookup.
ENSEMBL_SPECIES = {
    "Solanum lycopersicum": "solanum_lycopersicum",
    "Cucumis sativus": "cucumis_sativus",
    "Cucumis melo": "cucumis_melo",
    "Brassica oleracea": "brassica_oleracea",
    "Daucus carota": "daucus_carota",
    "Lactuca sativa": "lactuca_sativa",
}

RATE = float(os.environ.get("ATLAS_RATE_DELAY", "0.4"))
_last_call = [0.0]


# --------------------------------------------------------------------------
# HTTP with cache + rate limit
# --------------------------------------------------------------------------

def fetch(url: str, *, cache: bool = True) -> str | None:
    """GET a URL, caching the body. Returns None on any failure (never raises)."""
    CACHE.mkdir(exist_ok=True)
    key = hashlib.sha256(url.encode()).hexdigest()[:24]
    path = CACHE / f"{key}.txt"
    if cache and path.exists():
        return path.read_text(encoding="utf-8")

    wait = RATE - (time.time() - _last_call[0])
    if wait > 0:
        time.sleep(wait)
    _last_call[0] = time.time()

    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            body = r.read().decode("utf-8", "replace")
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as e:
        print(f"    ! network: {type(e).__name__}: {e}", file=sys.stderr)
        return None
    if cache:
        path.write_text(body, encoding="utf-8")
    return body


def fetch_json(url: str) -> dict | list | None:
    body = fetch(url)
    if not body:
        return None
    try:
        return json.loads(body)
    except json.JSONDecodeError:
        return None


def eutils_key() -> str:
    k = os.environ.get("NCBI_API_KEY", "")
    return f"&api_key={k}" if k else ""


# --------------------------------------------------------------------------
# Gene name handling
# --------------------------------------------------------------------------

def split_genes(field: str) -> list[str]:
    """'SlGAD3, SlGAD2' -> ['SlGAD3', 'SlGAD2'];  'Pun1 (AT3)' -> ['Pun1', 'AT3']."""
    if not field:
        return []
    field = re.sub(r"\b(plus|and|vs\.?|family|Network|candidates?)\b", ",", field, flags=re.I)
    parts = re.split(r"[,/;]| \+ ", field)
    out: list[str] = []
    for p in parts:
        p = p.strip().strip(".")
        if not p:
            continue
        # keep both 'Pun1' and the parenthesised synonym 'AT3'
        m = re.match(r"^(.*?)\s*\((.+?)\)\s*$", p)
        if m:
            for frag in (m.group(1), m.group(2)):
                frag = frag.strip()
                if frag and len(frag) <= 24 and not frag.lower().startswith("for "):
                    out.append(frag)
        elif len(p) <= 24:
            out.append(p)
    seen, uniq = set(), []
    for g in out:
        if g.lower() not in seen:
            seen.add(g.lower())
            uniq.append(g)
    return uniq


def parse_reference(ref: str) -> tuple[str, str]:
    """'Nonaka et al., 2017 (Sci. Rep.)' -> ('Nonaka', '2017'). Empty strings if absent."""
    if not ref:
        return "", ""
    author = ""
    m = re.match(r"\s*([A-Z][A-Za-zÀ-ɏ'\-]+)", ref)
    if m:
        author = m.group(1)
    year = ""
    y = re.search(r"\b(19[89]\d|20[0-4]\d)\b", ref)
    if y:
        year = y.group(1)
    return author, year


# --------------------------------------------------------------------------
# Database queries
# --------------------------------------------------------------------------

def ncbi_gene(symbol: str, species: str) -> dict:
    """Look a gene symbol up in NCBI Gene, constrained to one organism."""
    term = f'{symbol}[Gene Name] AND "{species}"[Organism]'
    url = (f"{EUTILS}/esearch.fcgi?db=gene&retmode=json&retmax=5"
           f"&term={urllib.parse.quote(term)}{eutils_key()}")
    res = fetch_json(url)
    ids = (res or {}).get("esearchresult", {}).get("idlist", []) if isinstance(res, dict) else []
    if not ids:
        return {"found": False, "query": term}
    url = f"{EUTILS}/esummary.fcgi?db=gene&retmode=json&id={','.join(ids)}{eutils_key()}"
    summ = fetch_json(url)
    hits = []
    if isinstance(summ, dict):
        for gid in ids:
            rec = summ.get("result", {}).get(gid)
            if not isinstance(rec, dict):
                continue
            hits.append({
                "ncbi_gene_id": gid,
                "official_symbol": rec.get("name", ""),
                "description": rec.get("description", ""),
                "other_aliases": rec.get("otheraliases", ""),
                "organism": (rec.get("organism") or {}).get("scientificname", ""),
                "genomic_context": (rec.get("genomicinfo") or [{}])[0].get("chrloc", ""),
            })
    return {"found": bool(hits), "query": term, "hits": hits}


def pubmed_reference(gene: str, author: str, year: str, trait: str) -> dict:
    """Find the PMID for a cited reference. Only returns what PubMed actually matched."""
    if not author and not gene:
        return {"found": False, "reason": "no searchable reference"}
    bits = []
    if author:
        bits.append(f"{author}[Author]")
    if year:
        bits.append(f"{year}[dp]")
    if gene:
        bits.append(f"{gene}[tiab]")
    term = " AND ".join(bits)
    url = (f"{EUTILS}/esearch.fcgi?db=pubmed&retmode=json&retmax=5"
           f"&term={urllib.parse.quote(term)}{eutils_key()}")
    res = fetch_json(url)
    ids = (res or {}).get("esearchresult", {}).get("idlist", []) if isinstance(res, dict) else []

    if not ids and gene and trait:
        # fall back to a topic search; flagged as weaker evidence of the exact citation
        term = f"{gene}[tiab] AND {trait}[tiab]"
        url = (f"{EUTILS}/esearch.fcgi?db=pubmed&retmode=json&retmax=5"
               f"&term={urllib.parse.quote(term)}{eutils_key()}")
        res = fetch_json(url)
        ids = (res or {}).get("esearchresult", {}).get("idlist", []) if isinstance(res, dict) else []
        if ids:
            return _pubmed_summaries(ids, term, exact=False)
        return {"found": False, "query": term}
    if not ids:
        return {"found": False, "query": term}
    return _pubmed_summaries(ids, term, exact=bool(author and year))


def _pubmed_summaries(ids: list[str], term: str, *, exact: bool) -> dict:
    url = f"{EUTILS}/esummary.fcgi?db=pubmed&retmode=json&id={','.join(ids)}{eutils_key()}"
    summ = fetch_json(url)
    hits = []
    if isinstance(summ, dict):
        for pid in ids:
            rec = summ.get("result", {}).get(pid)
            if not isinstance(rec, dict):
                continue
            doi = ""
            for aid in rec.get("articleids", []) or []:
                if aid.get("idtype") == "doi":
                    doi = aid.get("value", "")
            hits.append({
                "pmid": pid,
                "title": rec.get("title", ""),
                "journal": rec.get("source", ""),
                "year": (rec.get("pubdate", "") or "")[:4],
                "doi": doi,
            })
    return {
        "found": bool(hits),
        "query": term,
        "match_type": "author+year+gene" if exact else "gene+trait topic search",
        "hits": hits,
    }


def ensembl_gene(gene_id: str, symbol: str, species: str) -> dict:
    """Confirm a stable gene ID, or look the symbol up in Ensembl Plants."""
    if gene_id:
        cand = re.findall(r"[A-Za-z]{2,}[0-9][A-Za-z0-9._\-]+", gene_id)
        for c in cand:
            res = fetch_json(f"{ENSEMBL}/lookup/id/{urllib.parse.quote(c)}?content-type=application/json")
            if isinstance(res, dict) and res.get("id"):
                return {
                    "found": True, "source": "lookup/id", "ensembl_id": res.get("id"),
                    "biotype": res.get("biotype", ""), "species": res.get("species", ""),
                    "seq_region": res.get("seq_region_name", ""),
                    "start": res.get("start"), "end": res.get("end"),
                }
    key = ENSEMBL_SPECIES.get(species)
    if key and symbol:
        res = fetch_json(f"{ENSEMBL}/lookup/symbol/{key}/{urllib.parse.quote(symbol)}"
                         f"?content-type=application/json")
        if isinstance(res, dict) and res.get("id"):
            return {
                "found": True, "source": "lookup/symbol", "ensembl_id": res.get("id"),
                "biotype": res.get("biotype", ""), "species": res.get("species", ""),
                "seq_region": res.get("seq_region_name", ""),
                "start": res.get("start"), "end": res.get("end"),
            }
    return {"found": False}


# --------------------------------------------------------------------------
# Migration: v2 -> v3 skeleton (offline)
# --------------------------------------------------------------------------

CURATED_FIELDS = (
    "linkage_drag", "pleiotropic_effects", "tissue", "environmental_interaction",
    "conventional_breeding", "recommended_strategy", "editing_possibility",
    "trait_inheritance", "gene_action",
)


def migrate() -> dict:
    v2 = json.loads(V2.read_text(encoding="utf-8"))
    entries = []
    for e in v2.get("entries", []):
        genes = split_genes(e.get("gene", ""))
        entries.append({
            "id": e.get("id"),
            "crop": e.get("crop"),
            "scientific_name": e.get("scientific_name") or SPECIES.get(e.get("crop", ""), ""),
            "category": e.get("category"),
            "trait": e.get("trait"),
            "gene_label": e.get("gene"),
            "genes": [{"symbol": g, "claimed_gene_id": "", "verification": None} for g in genes],
            "portability": e.get("portability"),
            # claims carried over from v2, explicitly marked as not yet verified
            "claims": {
                "gene_function": e.get("gene_function", ""),
                "gene_id": e.get("gene_id", ""),
                "evidence_level": e.get("evidence_level", ""),
                "evidence_note": e.get("evidence_note", ""),
                "scientific_reference": e.get("scientific_reference", ""),
                "pmid_doi": e.get("pmid_doi", ""),
                "patent_reference": e.get("patent_reference", ""),
                **{k: e.get(k, "") for k in CURATED_FIELDS},
            },
            "verification": {
                "status": "unverified",
                "verified_on": None,
                "verified_by": None,
                "gene_ids": [],
                "pmids": [],
                "dois": [],
                "notes": [],
            },
        })
    out = {
        "_schema": {
            "version": "3.0.0",
            "platform": "Trait Transfer & Genome Editing Intelligence Platform",
            "principle": (
                "Gene identity and literature references are machine-verified against NCBI Gene, "
                "PubMed and Ensembl Plants. Portability, linkage drag, pleiotropy and IP remain "
                "indicative research signals curated by a breeder, not database facts."
            ),
            "verification_states": {
                "verified": "Gene ID and at least one reference confirmed in a public database.",
                "partial": "Gene confirmed, reference not resolved (or the reverse).",
                "unverified": "Not yet checked, or database returned no match.",
                "conflicted": "Database answer disagrees with the curated claim — needs human review.",
            },
            "evidence_levels": v2.get("_schema", {}).get("evidence_levels", {}),
            "portability_tiers": v2.get("_schema", {}).get("portability_tiers", {}),
        },
        "stats": {"total_entries": len(entries), "verified": 0, "partial": 0, "unverified": len(entries)},
        "entries": entries,
    }
    V3.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"migrated {len(entries)} entries -> {V3.relative_to(ROOT)}")
    return out


# --------------------------------------------------------------------------
# Verification pass
# --------------------------------------------------------------------------

def verify(only: str | None = None) -> dict:
    atlas = json.loads(V3.read_text(encoding="utf-8")) if V3.exists() else migrate()
    today = time.strftime("%Y-%m-%d")
    rows = [("entry_id", "crop", "gene", "ncbi_gene_id", "official_symbol",
             "ensembl_id", "pmid", "doi", "match_type", "status")]

    for e in atlas["entries"]:
        if only and only.lower() not in (e.get("crop") or "").lower():
            continue
        species = e.get("scientific_name") or ""
        author, year = parse_reference(e["claims"].get("scientific_reference", ""))
        trait = (e.get("trait") or "").split("(")[0].strip()
        print(f"  {e['id']}  {e.get('crop')}  [{len(e['genes'])} gene(s)]")

        gene_ids, pmids, dois, notes = [], [], [], []
        for g in e["genes"]:
            sym = g["symbol"]
            ncbi = ncbi_gene(sym, species) if species else {"found": False, "reason": "no species"}
            ens = ensembl_gene(e["claims"].get("gene_id", ""), sym, species)
            ref = pubmed_reference(sym, author, year, trait)

            g["verification"] = {"ncbi_gene": ncbi, "ensembl": ens, "pubmed": ref,
                                 "checked_on": today}

            ncbi_id = official = ens_id = pmid = doi = mtype = ""
            if ncbi.get("found"):
                h = ncbi["hits"][0]
                ncbi_id, official = h["ncbi_gene_id"], h["official_symbol"]
                gene_ids.append(f"NCBI:{ncbi_id}")
            if ens.get("found"):
                ens_id = ens["ensembl_id"]
                gene_ids.append(f"Ensembl:{ens_id}")
            if ref.get("found"):
                h = ref["hits"][0]
                pmid, doi, mtype = h["pmid"], h["doi"], ref.get("match_type", "")
                pmids.append(pmid)
                if doi:
                    dois.append(doi)

            gene_ok = bool(ncbi_id or ens_id)
            ref_ok = bool(pmid) and mtype == "author+year+gene"
            status = ("verified" if gene_ok and ref_ok else
                      "partial" if gene_ok or pmid else "unverified")
            if ncbi_id and official and official.lower() != sym.lower():
                notes.append(f"{sym}: NCBI official symbol is '{official}' — check synonym")
            if pmid and mtype != "author+year+gene":
                notes.append(f"{sym}: PMID {pmid} found by topic search, not by the cited "
                             f"author/year — citation not confirmed")
            rows.append((e["id"], e.get("crop", ""), sym, ncbi_id, official,
                         ens_id, pmid, doi, mtype, status))

        gene_ids = sorted(set(gene_ids))
        pmids = sorted(set(pmids))
        any_strong = any(
            (g.get("verification") or {}).get("pubmed", {}).get("match_type") == "author+year+gene"
            for g in e["genes"])
        status = ("verified" if gene_ids and pmids and any_strong else
                  "partial" if gene_ids or pmids else "unverified")
        e["verification"] = {
            "status": status,
            "verified_on": today,
            "verified_by": "verify_genes.py (NCBI Gene + PubMed + Ensembl Plants)",
            "gene_ids": gene_ids,
            "pmids": pmids,
            "dois": sorted(set(dois)),
            "notes": notes,
        }

    counts = {"verified": 0, "partial": 0, "unverified": 0, "conflicted": 0}
    for e in atlas["entries"]:
        counts[e["verification"]["status"]] = counts.get(e["verification"]["status"], 0) + 1
    atlas["stats"].update({"total_entries": len(atlas["entries"]), **counts,
                           "last_verified": today})
    V3.write_text(json.dumps(atlas, indent=2, ensure_ascii=False), encoding="utf-8")
    REPORT.write_text("\n".join("\t".join(str(c) for c in r) for r in rows) + "\n",
                      encoding="utf-8")
    print(f"\nwrote {V3.relative_to(ROOT)} and {REPORT.relative_to(ROOT)}")
    report(atlas)
    return atlas


def report(atlas: dict | None = None) -> None:
    atlas = atlas or json.loads(V3.read_text(encoding="utf-8"))
    s = atlas.get("stats", {})
    total = s.get("total_entries", 0)
    print("\n  status      entries")
    print("  " + "-" * 22)
    for k in ("verified", "partial", "unverified", "conflicted"):
        if k in s:
            print(f"  {k:<11} {s[k]:>5}")
    print(f"  {'total':<11} {total:>5}")
    flagged = [(e["id"], n) for e in atlas["entries"] for n in e["verification"].get("notes", [])]
    if flagged:
        print(f"\n  {len(flagged)} note(s) needing human review:")
        for eid, n in flagged[:40]:
            print(f"    {eid}: {n}")


def main() -> int:
    ap = argparse.ArgumentParser(description="Verify Transfer Atlas genes and references.")
    ap.add_argument("--migrate", action="store_true", help="build the v3 skeleton from v2 (offline)")
    ap.add_argument("--verify", action="store_true", help="query databases and fill verification")
    ap.add_argument("--report", action="store_true", help="print the current summary")
    ap.add_argument("--only", metavar="CROP", help="restrict to one crop")
    a = ap.parse_args()
    if not (a.migrate or a.verify or a.report):
        ap.print_help()
        return 1
    if a.migrate:
        migrate()
    if a.verify:
        verify(a.only)
    if a.report and not a.verify:
        report()
    return 0


if __name__ == "__main__":
    sys.exit(main())
