#!/usr/bin/env python3
"""Candidate-gene miner for the Transfer Atlas.

For a trait, pull PubMed abstracts that carry BOTH the trait and functional-
evidence language (overexpression, knockout, CRISPR, mutant...), then surface
the gene symbols that recur in them, flagged by causal phrasing. This is a
curation aid, not an answer: a human still decides which candidates are real
and fills the Atlas schema (portability, linkage, side-effects, IP).
"""
import collections
import csv
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
API_KEY = os.environ.get("NCBI_API_KEY", "")
DELAY = 0.12 if API_KEY else 0.4
MAX_ABS = int(os.environ.get("MAX_ABSTRACTS", "250"))
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results"

# Gene-symbol-like tokens: letters+digit core, optional species prefix / suffix (SlSOS1, OsHKT1;5, DREB2A)
GENE_RE = re.compile(r"\b([A-Za-z]{2,6}[0-9][A-Za-z0-9;.\-]{0,6})\b")
CAUSAL = re.compile(r"\b(improv|enhanc|confer|increas|reduc|require|regulat|trigger|restor|confers|mediat)", re.I)
STOP = {"CO2", "H2O", "H2O2", "NA+", "K+", "NACL", "ABA", "ROS", "ATP", "NADPH", "PCR", "QRT",
        "RNA", "DNA", "CDNA", "MRNA", "RT", "QPCR", "GWAS", "QTL", "SNP", "LED", "UV", "PH",
        "3D", "2D", "P5CS", "F1", "F2", "T2", "T3"}  # noisy non-genes (keep real ones like SOS1)


def call(endpoint, params, retries=4):
    params = dict(params, tool="transfer-atlas")
    if API_KEY:
        params["api_key"] = API_KEY
    url = BASE + endpoint + "?" + urllib.parse.urlencode(params)
    for attempt in range(retries):
        try:
            time.sleep(DELAY)
            with urllib.request.urlopen(url, timeout=90) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            time.sleep(2 ** attempt)
            last = e
    raise RuntimeError(f"Failed {url}: {last}")


def esearch(term, retmax):
    import json as _j
    d = _j.loads(call("esearch.fcgi", {"db": "pubmed", "term": term, "retmode": "json",
                                       "retmax": retmax, "sort": "relevance"}))
    return d["esearchresult"]["idlist"], int(d["esearchresult"]["count"])


def efetch(ids):
    import xml.etree.ElementTree as ET
    recs = []
    for i in range(0, len(ids), 100):
        xml = call("efetch.fcgi", {"db": "pubmed", "id": ",".join(ids[i:i + 100]), "retmode": "xml"})
        root = ET.fromstring(xml)
        for art in root.findall(".//PubmedArticle"):
            pmid = art.findtext(".//PMID") or ""
            title = "".join(art.find(".//ArticleTitle").itertext()) if art.find(".//ArticleTitle") is not None else ""
            abstract = " ".join("".join(t.itertext()) for t in art.findall(".//Abstract/AbstractText"))
            recs.append({"pmid": pmid, "title": title, "abstract": abstract})
    return recs


def main():
    trait = sys.argv[1] if len(sys.argv) > 1 else "salt"
    cfg = json.loads((ROOT / "traits.json").read_text())[trait]
    OUT.mkdir(exist_ok=True)
    term = f'{cfg["trait_terms"]} AND {cfg["evidence_terms"]} AND {cfg["context_terms"]}'
    ids, total = esearch(term, MAX_ABS)
    recs = efetch(ids)
    print(f"[{trait}] {total} functional-evidence hits; analysing {len(recs)} abstracts")

    gene_pmids = collections.defaultdict(set)
    gene_causal = collections.defaultdict(int)
    for r in recs:
        text = r["title"] + ". " + r["abstract"]
        causal_here = bool(CAUSAL.search(text))
        for m in set(GENE_RE.findall(text)):
            key = m.upper()
            if key in STOP or len(m) < 3:
                continue
            gene_pmids[m].add(r["pmid"])
            if causal_here:
                gene_causal[m] += 1

    rows = []
    for gene, pmids in gene_pmids.items():
        if len(pmids) < 2:  # need recurrence
            continue
        rows.append({"gene": gene, "papers": len(pmids), "causal_papers": gene_causal[gene],
                     "sample_pmids": ";".join(sorted(pmids)[:4])})
    rows.sort(key=lambda r: (-r["causal_papers"], -r["papers"]))

    out = OUT / f"candidates_{trait}.tsv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene", "papers", "causal_papers", "sample_pmids"], delimiter="\t")
        w.writeheader()
        w.writerows(rows[:120])
    print(f"[{trait}] wrote {min(len(rows),120)} candidate genes -> {out.name}")
    for r in rows[:25]:
        print(f"  {r['gene']:<12} papers={r['papers']:<3} causal={r['causal_papers']:<3} {r['sample_pmids']}")


if __name__ == "__main__":
    main()
