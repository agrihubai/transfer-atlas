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
  NCBI Gene      (E-utilities esearch/esummary)  -> gene ID, official symbol, description,
                                                    genomic coordinates, and cross-check of the
                                                    existence check of claimed locus IDs
  NCBI PubMed    (E-utilities esearch/esummary)  -> PMID + title + DOI for the cited reference

Ensembl Plants is not queried: its REST service lives on rest.ensemblgenomes.org, which this
environment cannot reach. rest.ensembl.org serves vertebrates only and rejects plant species.
Locus IDs are therefore confirmed against NCBI Gene instead.

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

# Species prefixes used in the literature but not in NCBI gene names (SlGAD3 -> GAD3).
PREFIXES = ("Sl", "Ca", "Cs", "Cm", "Cl", "Sm", "Ls", "Bo", "Br", "Dc", "Ac", "Le",
            "At", "Os", "Ta", "Hv", "Zm", "Gm", "Vv")

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
        g = g.strip(" ()[].")
        # a gene symbol starts with a letter and is more than one character;
        # this drops fragments left by splitting ('SOS1/2' -> '2', 'family)' -> ')')
        if len(g) < 2 or not re.match(r"^[A-Za-z]", g) or not re.search(r"[A-Za-z]{2}", g):
            continue
        if g.lower() in {"plus", "and", "family", "network", "candidate", "candidates",
                         "scavengers", "regulators", "catabolism", "related", "loci"}:
            continue
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

def strip_prefix(symbol: str) -> str:
    """SlGAD3 -> GAD3. Returns '' when the symbol carries no species prefix."""
    m = re.match(r"^([A-Z][a-z])([A-Z0-9].*)$", symbol)
    if m and m.group(1) in PREFIXES and len(m.group(2)) >= 2:
        return m.group(2)
    return ""


def _gene_search(term: str) -> list[str]:
    url = (f"{EUTILS}/esearch.fcgi?db=gene&retmode=json&retmax=10"
           f"&term={urllib.parse.quote(term)}{eutils_key()}")
    res = fetch_json(url)
    if not isinstance(res, dict):
        return []
    return res.get("esearchresult", {}).get("idlist", []) or []


def ncbi_gene(symbol: str, species: str, aliases: list[str] | None = None) -> dict:
    """Look a gene symbol up in NCBI Gene, constrained to one organism.

    Literature symbols carry a species prefix that NCBI does not use, so the
    prefix-stripped form is tried as well. Hits are ranked: a gene whose official
    symbol equals the query outranks one that merely lists it as an alias, which
    filters out unrelated genes that happen to share an alias.
    """
    queries, ids = [], []
    for sym in [symbol, strip_prefix(symbol), *(aliases or [])]:
        if not sym:
            continue
        term = f'{sym}[Gene Name] AND "{species}"[Organism]'
        queries.append(term)
        for i in _gene_search(term):
            if i not in ids:
                ids.append(i)
    if not ids:
        return {"found": False, "queries": queries}

    url = f"{EUTILS}/esummary.fcgi?db=gene&retmode=json&id={','.join(ids[:10])}{eutils_key()}"
    summ = fetch_json(url)
    hits = []
    wanted = ({symbol.lower(), strip_prefix(symbol).lower()}
              | {a.lower() for a in (aliases or [])}) - {""}
    if isinstance(summ, dict):
        for gid in ids[:10]:
            rec = summ.get("result", {}).get(gid)
            if not isinstance(rec, dict):
                continue
            name = rec.get("name", "") or ""
            aliases = [a.strip().lower() for a in (rec.get("otheraliases") or "").split(",")]
            gi = (rec.get("genomicinfo") or [{}])[0]
            if name.lower() in wanted:
                rank, how = 0, "official symbol"
            elif wanted & set(aliases):
                rank, how = 1, "listed as alias"
            else:
                rank, how = 2, "weak — symbol not in name or aliases"
            hits.append({
                "ncbi_gene_id": gid,
                "official_symbol": name,
                "description": rec.get("description", ""),
                "other_aliases": rec.get("otheraliases", ""),
                "organism": (rec.get("organism") or {}).get("scientificname", ""),
                "chromosome": gi.get("chrloc", "") or rec.get("chromosome", ""),
                "coordinates": (f"{gi.get('chraccver','')}:{gi.get('chrstart','')}-{gi.get('chrstop','')}"
                                if gi.get("chraccver") else ""),
                "match": how,
                "_rank": rank,
            })
    hits.sort(key=lambda h: h["_rank"])
    # drop weak hits when a strong one exists, so the first hit is never a coincidence
    if hits and hits[0]["_rank"] <= 1:
        hits = [h for h in hits if h["_rank"] <= 1]
    for h in hits:
        h.pop("_rank", None)
    return {"found": bool(hits), "queries": queries, "hits": hits}


def locus_genbank(locus: str, species: str) -> dict:
    """For a locus tag that NCBI holds in a sequence record, read what the record says.

    Plant genes absent from NCBI Gene are still identifiable this way: a GenBank
    record whose /locus_tag is the tag, carrying /product and /protein_id, anchors
    the identity as well as a Gene record does.
    """
    term = f'{locus} AND "{species}"[Organism]' if species else locus
    url = (f"{EUTILS}/esearch.fcgi?db=nuccore&retmode=json&retmax=2"
           f"&term={urllib.parse.quote(term)}{eutils_key()}")
    res = fetch_json(url)
    ids = ((res or {}).get("esearchresult", {}) or {}).get("idlist", []) if isinstance(res, dict) else []
    if not ids:
        return {}
    gb = fetch(f"{EUTILS}/efetch.fcgi?db=nuccore&id={ids[0]}&rettype=gb&retmode=text"
               f"{eutils_key()}")
    if not gb or locus not in gb:
        return {}
    def grab(key):
        m = re.search(r'/' + key + r'="([^"]+)"', gb)
        return m.group(1) if m else ""
    acc = re.search(r"^ACCESSION\s+(\S+)", gb, re.M)
    if grab("locus_tag") != locus:
        return {}
    return {"accession": acc.group(1) if acc else "", "product": grab("product"),
            "protein_id": grab("protein_id")}


def genbank_accession(acc: str, species: str) -> dict:
    """Confirm a GenBank accession pinned to a gene: organism, gene name, product.

    Used where neither NCBI Gene nor a locus tag resolves the symbol. A wrong
    accession fails here instead of being trusted: the organism must match.
    """
    if not acc:
        return {}
    summ = fetch_json(f"{EUTILS}/esummary.fcgi?db=nuccore&retmode=json"
                      f"&id={urllib.parse.quote(acc)}{eutils_key()}")
    rec = None
    if isinstance(summ, dict):
        for k in (summ.get("result", {}) or {}).get("uids", []) or []:
            rec = summ["result"].get(k)
    if not isinstance(rec, dict) or not rec.get("title"):
        return {"error": f"accession {acc} did not resolve"}
    org = (rec.get("organism") or "")
    title = rec.get("title", "")
    if species and species.split()[0].lower() not in (org + " " + title).lower():
        return {"error": f"accession {acc} is {org or title[:40]}, not {species}"}
    gb = fetch(f"{EUTILS}/efetch.fcgi?db=nuccore&id={urllib.parse.quote(acc)}"
               f"&rettype=gb&retmode=text{eutils_key()}") or ""
    def grab(key):
        m = re.search(r'/' + key + r'="([^"]+)"', gb)
        return m.group(1) if m else ""
    return {"accession": rec.get("accessionversion", acc), "title": title,
            "gene": grab("gene"), "product": grab("product"),
            "protein_id": grab("protein_id")}


def check_locus_id(locus_field: str, species: str) -> dict:
    """Try to confirm a claimed locus ID (Solyc01g005000, Csa5G623470, …).

    LIMITATION, stated rather than hidden. Genome-project locus tags — ITAG
    (Solyc), Csa, Cla97, MELO3C, SMEL4.1, BraA — are not indexed in NCBI at all:
    NCBI Gene carries its own symbols (MLO1) and numeric GeneIDs instead, and
    nucleotide->gene elinks do not resolve these tags either. So a 'no match' here
    says nothing about whether the ID is right; it says NCBI is the wrong place to
    ask. An UNCONSTRAINED text search does return hits, but they are incidental
    matches in unrelated organisms (searching Solyc04g049090 returns a Genlisea
    aurea scaffold), so the search is constrained by organism and a cross-species
    hit is never counted as confirmation.

    Resolving these IDs needs Ensembl Plants (rest.ensemblgenomes.org) or Sol
    Genomics Network (solgenomics.net); neither is reachable from this
    environment. Until one is, the locus-to-gene mapping stays a curated claim
    and is labelled as such.
    """
    loci = re.findall(r"[A-Za-z]{3,}\d[A-Za-z0-9._\-]{4,}", locus_field or "")[:6]
    if not loci:
        return {"checked": False}
    results = []
    for locus in loci:
        term = f'{locus} AND "{species}"[Organism]' if species else locus
        url = (f"{EUTILS}/esearch.fcgi?db=nuccore&retmode=json&retmax=3"
               f"&term={urllib.parse.quote(term)}{eutils_key()}")
        res = fetch_json(url)
        sr = (res or {}).get("esearchresult", {}) if isinstance(res, dict) else {}
        count = int(sr.get("count", "0") or 0)
        gb = locus_genbank(locus, species) if count else {}
        results.append({
            "locus": locus,
            "query": term,
            "genbank_identity": gb or None,
            "ncbi_records_in_species": count,
            "state": "confirmed_in_species" if count else "not_checkable_here",
            "status": (f"appears in {count} NCBI sequence record(s) for {species}"
                       if count else
                       "not indexed in NCBI — needs Ensembl Plants or Sol Genomics; "
                       "this is not evidence the ID is wrong"),
        })
    return {
        "checked": True,
        "scope": "identifier presence within the stated species; never gene assignment",
        "results": results,
    }


def pubmed_reference(gene: str, author: str, year: str, trait: str,
                     crop: str = "", journal: str = "") -> dict:
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

    if not ids and author and year:
        # The cited author and year are often right while the gene symbol never
        # appears in the title or abstract as the draft spells it. Search on
        # author+year and pick the paper whose TITLE carries the crop and the
        # trait. A surname can be shared across fields entirely — "Nekrasov 2017"
        # returns cardiology and oncology papers — so a candidate is accepted
        # only on term overlap, never because it came first.
        cand = _author_year_titles(author, year)
        best = _pick_by_title(cand, trait, crop, gene, journal)
        if best and best.get("near_miss"):
            return {"found": True, "query": f"{author}[Author] AND {year}[dp]",
                    "match_type": "awaiting human confirmation: journal and author "
                                  "match, title does not confirm the trait",
                    "matched_terms": best["terms"], "hits": [best["hit"]]}
        if best:
            return {"found": True, "query": f"{author}[Author] AND {year}[dp]",
                    "match_type": "author+year+journal, trait matched in title",
                    "matched_terms": best["terms"], "candidates_considered": len(cand),
                    "hits": [best["hit"]]}

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


STOP = {"and", "the", "for", "with", "via", "under", "from", "into", "high", "low",
        "fruit", "plant", "gene", "trait", "type", "development", "accumulation"}


def _author_year_titles(author: str, year: str) -> list[dict]:
    term = f"{author}[Author] AND {year}[dp]"
    url = (f"{EUTILS}/esearch.fcgi?db=pubmed&retmode=json&retmax=20"
           f"&term={urllib.parse.quote(term)}{eutils_key()}")
    res = fetch_json(url)
    ids = ((res or {}).get("esearchresult", {}) or {}).get("idlist", []) if isinstance(res, dict) else []
    if not ids:
        return []
    summ = fetch_json(f"{EUTILS}/esummary.fcgi?db=pubmed&retmode=json"
                      f"&id={','.join(ids)}{eutils_key()}")
    out = []
    if isinstance(summ, dict):
        for pid in ids:
            rec = summ.get("result", {}).get(pid)
            if not isinstance(rec, dict):
                continue
            doi = ""
            for aid in rec.get("articleids", []) or []:
                if aid.get("idtype") == "doi":
                    doi = aid.get("value", "")
            out.append({"pmid": pid, "title": rec.get("title", ""),
                        "journal": rec.get("source", ""),
                        "year": (rec.get("pubdate", "") or "")[:4], "doi": doi})
    return out


def parse_journal(ref: str) -> str:
    """Pull the journal out of a draft citation: 'Uluisik et al., 2016, Nat. Biotechnol.'"""
    if not ref:
        return ""
    m = re.search(r"\b(19|20)\d{2}\b[,.)\s]*\(?([A-Z][A-Za-z.&\s]{3,40})", ref)
    if not m:
        return ""
    return re.sub(r"[.\s]+$", "", m.group(2)).strip()


def _journal_matches(cited: str, candidate: str) -> bool:
    """Do two journal names refer to the same journal, allowing for abbreviation?"""
    def toks(x):
        x = (x or "").lower().replace(".", " ")
        return [t for t in re.findall(r"[a-z]+", x) if t not in {"of", "the", "and", "in", "for"}]
    a, b = toks(cited), toks(candidate)
    if not a or not b:
        return False
    # every cited token must prefix-match a candidate token: 'nat biotechnol' ~ 'Nat Biotechnol',
    # 'new phytol' ~ 'New Phytologist'
    used = set()
    for t in a:
        hit = next((i for i, u in enumerate(b)
                    if i not in used and (u.startswith(t) or t.startswith(u))), None)
        if hit is None:
            return False
        used.add(hit)
    return True


def _pick_by_title(cands: list[dict], trait: str, crop: str, gene: str,
                   cited_journal: str) -> dict | None:
    """Choose the candidate that is the cited paper, or none at all.

    Three guards, each added after a wrong match got through:
    - the journal named in the draft citation must match the candidate's journal.
      Without this, 'Uluisik 2016' matched an immunohistochemistry paper and
      'Traka 2013' matched a human nutrition trial by the same author in the same
      year as the genetics paper actually cited;
    - gene tokens under three characters are ignored, since 'PL' matches inside
      unrelated words;
    - every candidate is scored and the best is taken, not the first.
    A citation that names no journal is left unconfirmed: author and year alone
    cannot tell two papers apart, and guessing is worse than an open gap.
    """
    if not cited_journal:
        return None
    traitw = {w for w in re.findall(r"[a-z]{4,}", (trait or "").lower()) if w not in STOP}
    cropw = {w for w in re.findall(r"[a-z]{4,}", (crop or "").lower()) if w not in STOP}
    genew = {w for w in ({(gene or "").lower(), strip_prefix(gene or "").lower()} - {""})
             if len(w) >= 3}
    best = None
    near = None
    for c in cands:
        if not _journal_matches(cited_journal, c.get("journal", "")):
            continue
        t = (c["title"] or "").lower()
        # a corrigendum or erratum carries the same title as the paper it fixes;
        # prefer the paper itself
        if re.match(r"\s*(corrigendum|erratum|retraction|author correction)", t):
            if near is None:
                near = {"hit": c, "terms": [], "deprioritised": "corrigendum/erratum"}
            continue
        def found(words):
            return {w for w in words if re.search(r"\b" + re.escape(w), t)}
        hit_t, hit_c, hit_g = found(traitw), found(cropw), found(genew)
        score = 2 * len(hit_g) + len(hit_c) + len(hit_t)
        if len(hit_g | hit_c | hit_t) < 2:
            # journal agrees but the title does not carry enough of the trait to
            # confirm. Not evidence — but worth showing a human, who can settle
            # it in one look. See near_miss below.
            if near is None or near.get("deprioritised"):
                near = {"hit": c, "terms": sorted(hit_g | hit_c | hit_t)}
            continue
        if best is None or score > best["score"]:
            best = {"hit": c, "terms": sorted(hit_g | hit_c | hit_t), "score": score}
    if best is None and near is not None and not near.get("deprioritised"):
        near["near_miss"] = True
    if best is None and near is not None:
        near["near_miss"] = True
        return near
    return best


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


# --------------------------------------------------------------------------
# Migration: v2 -> v3 skeleton (offline)
# --------------------------------------------------------------------------

CURATED_FIELDS = (
    "linkage_drag", "pleiotropic_effects", "tissue", "environmental_interaction",
    "conventional_breeding", "recommended_strategy", "editing_possibility",
    "trait_inheritance", "gene_action",
)


CORRECTIONS = DATA / "corrections.json"


def load_corrections() -> dict:
    if not CORRECTIONS.exists():
        return {}
    try:
        return json.loads(CORRECTIONS.read_text(encoding="utf-8")).get("entries", {}) or {}
    except json.JSONDecodeError:
        return {}


def migrate() -> dict:
    """Build the working atlas from the untouched draft plus curated corrections.

    The draft (atlas_v2.json) is never edited, so the difference between what an
    external tool proposed and what a breeder decided stays visible in git and in
    the corrections file.
    """
    v2 = json.loads(V2.read_text(encoding="utf-8"))
    fixes = load_corrections()
    entries = []
    for e in v2.get("entries", []):
        fix = fixes.get(e.get("id")) or {}
        if fix.get("set"):
            e = {**e, **fix["set"]}
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
                "commercial_product": e.get("commercial_product", ""),
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
            "correction": {k: v for k, v in fix.items() if k != "set"} or None,
        })
    out = {
        "_schema": {
            "version": "3.0.0",
            "platform": "Trait Transfer & Genome Editing Intelligence Platform",
            "principle": (
                "Gene identity and literature references are machine-verified against NCBI Gene "
                "and PubMed. Portability, linkage drag, pleiotropy and IP remain "
                "indicative research signals curated by a breeder, not database facts."
            ),
            "verification_states": {
                "verified": "Gene ID and at least one reference confirmed in a public database.",
                "partial": "Gene confirmed, reference not resolved (or the reverse).",
                "unverified": "Not yet checked, or database returned no match.",
                "conflicted": "Database answer disagrees with the curated claim — needs human review.",
            },
            "evidence_levels": v2.get("_schema", {}).get("evidence_levels", {}),
            "level_derivation": DERIVATION_RULE,
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
             "locus_ids_exist", "pmid", "doi", "match_type", "status")]

    for e in atlas["entries"]:
        if only and only.lower() not in (e.get("crop") or "").lower():
            continue
        species = e.get("scientific_name") or ""
        author, year = parse_reference(e["claims"].get("scientific_reference", ""))
        trait = (e.get("trait") or "").split("(")[0].strip()
        print(f"  {e['id']}  {e.get('crop')}  [{len(e['genes'])} gene(s)]")

        gene_ids, pmids, dois, notes = [], [], [], []

        # Locus IDs are listed per entry, not per gene, so check them once.
        locus = check_locus_id(e["claims"].get("gene_id", ""), species)
        loci_exist = [r["locus"] for r in locus.get("results", [])
                      if r["state"] == "confirmed_in_species"]
        unresolved = [r["locus"] for r in locus.get("results", [])
                      if r["state"] == "not_checkable_here"]
        if unresolved:
            # a limitation of this pipeline, not a fault in the data
            notes.append(f"locus ID(s) {', '.join(unresolved)} not checkable against NCBI "
                         f"(genome-project tags are not indexed there) — needs Ensembl Plants "
                         f"or Sol Genomics; the IDs may well be correct")
        e["locus_check"] = locus

        for g in e["genes"]:
            sym = g["symbol"]
            corr = e.get("correction") or {}
            ncbi_id = official = pmid = doi = mtype = ""
            al = corr.get("search_aliases", {}).get(sym, [])
            ncbi = (ncbi_gene(sym, species, al) if species
                    else {"found": False, "reason": "no species"})
            # Identity can also come from a GenBank record carrying the locus tag,
            # which is how plant genes absent from NCBI Gene are pinned down.
            gb_id = {}
            mapped = corr.get("gene_locus", {}).get(sym, "")
            if not ncbi.get("found") and mapped:
                gb_id = locus_genbank(mapped, species)
                if gb_id:
                    gb_id["locus"] = mapped
            gid = corr.get("gene_ncbi_id", {}).get(sym, "")
            if not ncbi.get("found") and gid:
                summ = fetch_json(f"{EUTILS}/esummary.fcgi?db=gene&retmode=json"
                                  f"&id={gid}{eutils_key()}")
                rec = ((summ or {}).get("result", {}) or {}).get(str(gid)) \
                    if isinstance(summ, dict) else None
                org = (rec or {}).get("organism", {}).get("scientificname", "")
                if rec and (not species or org == species):
                    ncbi = {"found": True, "queries": [f"pinned id {gid}"], "hits": [{
                        "ncbi_gene_id": str(gid), "official_symbol": rec.get("name", ""),
                        "description": rec.get("description", ""),
                        "other_aliases": rec.get("otheraliases", ""), "organism": org,
                        "chromosome": rec.get("chromosome", ""), "coordinates": "",
                        "match": "pinned by curator"}]}
                elif rec:
                    notes.append(f"{sym}: pinned gene id {gid} is {org}, not {species}")
                else:
                    notes.append(f"{sym}: pinned gene id {gid} did not resolve")
            acc = corr.get("gene_accession", {}).get(sym, "")
            if not ncbi.get("found") and not gb_id and acc:
                got = genbank_accession(acc, species)
                if got.get("error"):
                    notes.append(f"{sym}: {got['error']}")
                elif got:
                    gb_id = got
            if gb_id and not gb_id.get("error"):
                gene_ids.append(f"GenBank:{gb_id.get('accession','')}")
            if ncbi.get("found"):
                h = ncbi["hits"][0]
                ncbi_id, official = h["ncbi_gene_id"], h["official_symbol"]
                gene_ids.append(f"NCBI:{ncbi_id}")
            ref = pubmed_reference(sym, author, year, trait,
                                   f"{e.get('crop','')} {species}",
                                   parse_journal(e["claims"].get("scientific_reference", "")))

            g["verification"] = {"ncbi_gene": ncbi, "genbank_locus": gb_id or None,
                                 "pubmed": ref, "checked_on": today}

            if ref.get("found"):
                h = ref["hits"][0]
                pmid, doi, mtype = h["pmid"], h["doi"], ref.get("match_type", "")
                pmids.append(pmid)
                if doi:
                    dois.append(doi)

            gene_ok = bool(ncbi_id or gb_id)
            ref_ok = mtype in ("author+year+gene",
                               "author+year+journal, trait matched in title")
            status = ("verified" if gene_ok and ref_ok else
                      "partial" if gene_ok or pmid else "unverified")
            if ncbi_id and ncbi["hits"][0]["match"] != "official symbol":
                notes.append(f"{sym}: matched NCBI gene {official} ({ncbi_id}) "
                             f"as {ncbi['hits'][0]['match']} — confirm it is the same gene")
            if pmid and mtype not in ("author+year+gene",
                                      "author+year+journal, trait matched in title"):
                notes.append(f"{sym}: PMID {pmid} found by topic search, not by the cited "
                             f"author/year — citation not confirmed")
            rows.append((e["id"], e.get("crop", ""), sym, ncbi_id, official,
                         ";".join(loci_exist), pmid, doi, mtype, status))

        gene_ids = sorted(set(gene_ids))
        pmids = sorted(set(pmids))
        any_strong = any(
            (g.get("verification") or {}).get("pubmed", {}).get("match_type")
            in ("author+year+gene", "author+year+journal, trait matched in title")
            for g in e["genes"])
        curated = check_curated(e["id"])
        curated_ok = [c for c in curated if c.get("pipeline_check", "").startswith("confirmed")
                      or c.get("pipeline_check", "").startswith("PMID resolves")]
        for c in curated_ok:
            if c.get("pmid") and c["pmid"] not in pmids:
                pmids.append(c["pmid"])
            if c.get("pubmed_doi") and c["pubmed_doi"] not in dois:
                dois.append(c["pubmed_doi"])

        n_genes = len(e["genes"])
        n_gene_ok = sum(1 for g in e["genes"]
                        if ((g.get("verification") or {}).get("ncbi_gene", {}).get("found")
                            or (g.get("verification") or {}).get("genbank_locus")))
        # "verified" requires every named gene to resolve, not just one of them —
        # otherwise an entry listing five genes passes on the strength of one.
        # A curated reference confirmed in PubMed counts as a real reference, but it
        # cannot by itself make an entry verified: every named gene must resolve too.
        ref_ok = any_strong or bool(curated_ok)
        status = ("verified" if n_genes and n_gene_ok == n_genes and ref_ok else
                  "partial" if gene_ids or pmids else "unverified")
        e["verification"] = {
            "status": status,
            "genes_resolved": f"{n_gene_ok}/{n_genes}",
            "curated_references": curated,
            "verified_on": today,
            "verified_by": "verify_genes.py (NCBI Gene + PubMed)",
            "gene_ids": gene_ids,
            "pmids": pmids,
            "dois": sorted(set(dois)),
            "notes": notes,
        }

    for e in atlas["entries"]:
        e["evidence"] = derive_level(e)

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


CURATED = DATA / "curated_references.json"


def check_curated(entry_id: str) -> list[dict]:
    """Check human- or tool-supplied references for one entry against PubMed.

    A suggestion is never taken on trust: the PMID must resolve, and the title,
    journal, year and DOI are taken from PubMed rather than from the suggestion.
    A DOI that disagrees with PubMed's is reported as a mismatch.
    """
    if not CURATED.exists():
        return []
    try:
        book = json.loads(CURATED.read_text(encoding="utf-8")).get("entries", {})
    except json.JSONDecodeError:
        return []
    out = []
    for ref in book.get(entry_id, []) or []:
        pmid = str(ref.get("pmid", "")).strip()
        rec = {k: ref.get(k, "") for k in
               ("pmid", "doi", "suggested_by", "supports", "does_not_support")}
        if not pmid:
            rec["pipeline_check"] = "no PMID given — not shown"
            out.append(rec)
            continue
        summ = fetch_json(f"{EUTILS}/esummary.fcgi?db=pubmed&retmode=json"
                          f"&id={pmid}{eutils_key()}")
        got = (summ or {}).get("result", {}).get(pmid) if isinstance(summ, dict) else None
        if not isinstance(got, dict) or not got.get("title"):
            rec["pipeline_check"] = "PMID did not resolve in PubMed — not shown"
            out.append(rec)
            continue
        doi = ""
        for aid in got.get("articleids", []) or []:
            if aid.get("idtype") == "doi":
                doi = aid.get("value", "")
        rec.update({
            "pipeline_check": "confirmed in PubMed",
            "pubmed_title": got.get("title", ""),
            "pubmed_journal": got.get("source", ""),
            "pubmed_year": (got.get("pubdate", "") or "")[:4],
            "pubmed_doi": doi,
        })
        if ref.get("doi") and doi and ref["doi"].lower() != doi.lower():
            rec["pipeline_check"] = (f"PMID resolves, but the suggested DOI does not match "
                                     f"PubMed's ({doi}) — suggestion corrected")
        out.append(rec)
    return out


def fetch_accessions() -> dict:
    """Attach a RefSeq transcript accession to every resolved gene.

    A BLAST deep link needs a query sequence. Without an accession the link only
    opens an empty form, so genes without one get no BLAST link in the UI rather
    than a link that does nothing.
    """
    atlas = json.loads(V3.read_text(encoding="utf-8"))
    n = 0
    for e in atlas["entries"]:
        for g in e["genes"]:
            v = g.get("verification") or {}
            hits = (v.get("ncbi_gene") or {}).get("hits") or []
            if not hits:
                continue
            uid = hits[0]["ncbi_gene_id"]
            acc = ""
            for linkname in ("gene_nuccore_refseqrna", "gene_nuccore_refseqgene"):
                url = (f"{EUTILS}/elink.fcgi?dbfrom=gene&db=nuccore&retmode=json"
                       f"&linkname={linkname}&id={uid}{eutils_key()}")
                res = fetch_json(url)
                ids = []
                if isinstance(res, dict):
                    for ls in res.get("linksets", []) or []:
                        for ldb in ls.get("linksetdbs", []) or []:
                            ids += ldb.get("links", []) or []
                if not ids:
                    continue
                summ = fetch_json(f"{EUTILS}/esummary.fcgi?db=nuccore&retmode=json"
                                  f"&id={ids[0]}{eutils_key()}")
                if isinstance(summ, dict):
                    rec = summ.get("result", {}).get(str(ids[0]), {})
                    acc = rec.get("accessionversion") or rec.get("caption") or ""
                if acc:
                    v["refseq_accession"] = acc
                    v["refseq_link_type"] = linkname
                    n += 1
                    break
    atlas["stats"]["genes_with_accession"] = n
    V3.write_text(json.dumps(atlas, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"attached {n} RefSeq accession(s)")
    return atlas


DERIVATION_RULE = (
    "A reference counts as confirmed only when PubMed matched the cited author, year and gene, "
    "or when a supplied reference resolved by its PMID. A PMID found by topic search alone does "
    "NOT count: it may be a different paper about the same gene. "
    "Level 1 is never derived — it requires a human to name a commercial cultivar or released "
    "product in corrections.json (commercial_product). "
    "Level 2 — every gene named in the entry resolves in the stated species AND a reference is "
    "confirmed. "
    "Level 3 — a reference is confirmed but gene identity is incomplete: either only some of the "
    "named genes resolve, or none of them do. "
    "Level 4 — genes resolve but no reference is confirmed: a candidate by identity only. "
    "Level 5 — neither a gene record nor a confirmed reference."
)


def derive_level(entry: dict) -> dict:
    """Compute an evidence level from what was verified, instead of trusting the draft."""
    genes = entry.get("genes") or []
    n = len(genes)
    ok = sum(1 for g in genes
             if (((g.get("verification") or {}).get("ncbi_gene") or {}).get("found")
                 or (g.get("verification") or {}).get("genbank_locus")))
    ver = entry.get("verification") or {}
    # Only confirmed references count. A topic-search PMID is not evidence for the
    # cited claim, so it must not raise the level.
    CONFIRMED = {"author+year+gene", "author+year+journal, trait matched in title"}
    strong = any(((g.get("verification") or {}).get("pubmed") or {}).get("match_type")
                 in CONFIRMED for g in genes)
    supplied = any(str(c.get("pipeline_check", "")).startswith(("confirmed", "PMID resolves"))
                   for c in (ver.get("curated_references") or []))
    ref = strong or supplied
    product = (entry.get("claims") or {}).get("commercial_product", "")

    if product:
        lvl, why = "Level 1", f"commercial product named by a human: {product}"
    elif n and ok == n and ref:
        lvl, why = "Level 2", f"all {n} gene(s) resolved and a reference confirmed"
    elif ref and ok:
        lvl, why = "Level 3", f"only {ok}/{n} genes resolved, with a confirmed reference"
    elif ref:
        lvl, why = ("Level 3", "a reference is confirmed, but none of the named gene symbols "
                               "resolve in the stated species — identity unestablished")
    elif ok:
        lvl, why = "Level 4", f"{ok}/{n} genes resolved, no confirmed reference"
    else:
        lvl, why = "Level 5", "no gene record and no confirmed reference resolved"

    claimed = (entry.get("claims") or {}).get("evidence_level", "")
    out = {"derived_level": lvl, "basis": why, "rule": "see _schema.level_derivation",
           "draft_claimed": claimed}
    if claimed and claimed != lvl:
        out["disagrees_with_draft"] = (
            f"the draft claimed {claimed}; this pipeline derives {lvl} from what it could verify")
    return out


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
    ap.add_argument("--accessions", action="store_true",
                    help="attach RefSeq accessions to resolved genes (for BLAST links)")
    ap.add_argument("--report", action="store_true", help="print the current summary")
    ap.add_argument("--only", metavar="CROP", help="restrict to one crop")
    a = ap.parse_args()
    if not (a.migrate or a.verify or a.accessions or a.report):
        ap.print_help()
        return 1
    if a.migrate:
        migrate()
    if a.verify:
        verify(a.only)
    if a.accessions:
        fetch_accessions()
    if a.report and not a.verify:
        report()
    return 0


if __name__ == "__main__":
    sys.exit(main())
