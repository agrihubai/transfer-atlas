#!/usr/bin/env python3
"""
Find a gene's likely counterpart in a target crop, by protein similarity.

Why this exists
---------------
Every transfer entry eventually asks the same question: the solution is
documented in species A, so which gene is it in crop B? Without an answer the
entry cannot be edited, and no sequence-based patent search is even possible.

How it answers it, and what the answer is worth
-----------------------------------------------
It aligns a set of query proteins against the target crop's whole proteome and
reports the best matches. The important part is that the queries include the
gene's CLOSE RELATIVES, not just the gene itself: a best match to one protein
means little when a paralogue would score almost as well. By aligning the
family together, each target protein can be labelled with which query it is
closest to — so a hit that is nearer the C-class paralogue than the D-class gene
of interest is visible as such instead of being reported as the ortholog.

This produces a ranked, reproducible candidate list. It does NOT establish
orthology: that needs a phylogenetic tree with the conserved domains aligned,
and ideally expression in the right tissue. The output says so, and the atlas
records these as candidates.

Usage
-----
  python3 scripts/find_ortholog.py --preset agl11-watermelon
  python3 scripts/find_ortholog.py --queries NP_001268118.1=VviAGL11,NP_850377.1=AtSHP2 \\
                                   --proteome <url> --label "Citrullus lanatus" --top 8
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
CACHE = ROOT / ".cache"
OUT = DATA / "ortholog_candidates.json"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
QUERY_TITLES: dict[str, str] = {}
PEP = "https://ftp.ebi.ac.uk/ensemblgenomes/pub/current/plants/fasta"
UA = "TransferAtlas/3.0 (https://agrihubai.github.io/transfer-atlas)"

PRESETS = {
    # The D-class gene of interest, its grape counterpart, and the C-class
    # paralogues it is most easily confused with.
    "agl11-watermelon": {
        "label": "Citrullus lanatus",
        "proteome": f"{PEP}/citrullus_lanatus/pep/Citrullus_lanatus.Cla97_v1.pep.all.fa.gz",
        "queries": {
            "NP_001268118.1": "VviAGL11 (grape, D-class, the entry's gene)",
            "NP_001329612.1": "AtSTK/AGL11 (Arabidopsis, D-class)",
            "NP_001325709.1": "AtSHP1/AGL1 (Arabidopsis, C-class paralogue)",
            "NP_850377.1": "AtSHP2/AGL5 (Arabidopsis, C-class paralogue)",
            "NP_001328877.1": "AtAG (Arabidopsis, C-class)",
        },
        "of_interest": ("VviAGL11", "AtSTK"),
        "entry": "WAT-SEED-01",
    },
}


_last = [0.0]


def fetch(url: str, binary: bool = False):
    """Cached GET. NCBI allows 3 requests a second and answers 429 when pushed
    harder, which is how an earlier version lost a whole job mid-run."""
    CACHE.mkdir(exist_ok=True)
    key = urllib.parse.quote(url, safe="")[-120:]
    path = CACHE / f"seq_{abs(hash(key)) % (10 ** 12)}{'.bin' if binary else '.txt'}"
    if path.exists():
        return path.read_bytes() if binary else path.read_text(encoding="utf-8")
    if "ncbi.nlm.nih.gov" in url:
        wait = 0.4 - (time.time() - _last[0])
        if wait > 0:
            time.sleep(wait)
        _last[0] = time.time()
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    last_err = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                raw = r.read()
            break
        except Exception as exc:   # 429 and reset connections are worth retrying
            last_err = exc
            time.sleep(2 + 4 * attempt)
    else:
        raise RuntimeError(f"giving up on {url[:70]}: {last_err}")
    if binary:
        path.write_bytes(raw)
        return raw
    text = raw.decode("utf-8", "replace")
    path.write_text(text, encoding="utf-8")
    return text


def query_sequences(accs: dict[str, str]) -> dict[str, tuple[str, str]]:
    """Fetch one query protein per id, keeping the real accession from the header.

    Fetched one at a time on purpose: ids may be numeric UIDs while the returned
    FASTA header carries an accession, so a batch fetch cannot be mapped back
    reliably — an earlier version matched id against header, silently got nothing,
    and reported an empty result as if the search had run.
    """
    got: dict[str, tuple[str, str]] = {}
    for ident, label in accs.items():
        fa = fetch(f"{EUTILS}/efetch.fcgi?db=protein&id={urllib.parse.quote(str(ident))}"
                   f"&rettype=fasta&retmode=text")
        lines = [l for l in (fa or "").splitlines() if l.strip()]
        if not lines or not lines[0].startswith(">"):
            print(f"  ! no sequence for {ident} ({label})", file=sys.stderr)
            continue
        acc = lines[0][1:].split()[0]
        title = " ".join(lines[0][1:].split()[1:])[:70]
        seq = "".join(l.strip() for l in lines[1:] if not l.startswith(">"))
        if not seq:
            print(f"  ! empty sequence for {ident} ({label})", file=sys.stderr)
            continue
        got[label] = (acc, seq)
        QUERY_TITLES[label] = f"{acc} · {len(seq)} aa · {title}"
        print(f"   {acc:<16}{len(seq):>5} aa  {label}")
        print(f"                    {title}")
    if not got:
        print("  ! no query sequences at all; refusing to report an empty search",
              file=sys.stderr)
    return got


def proteome(url: str) -> dict[str, str]:
    raw = fetch(url, binary=True)
    # Ensembl ships gzip; Sol Genomics ships plain FASTA.
    if url.endswith(".gz"):
        text = gzip.GzipFile(fileobj=io.BytesIO(raw)).read().decode("utf-8", "replace")
    else:
        text = raw.decode("utf-8", "replace")
    seqs, name, buf = {}, None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            # >ID pep ... gene:GENEID ...
            head = line[1:]
            first = head.split()[0]
            gid = next((p.split(":", 1)[1] for p in head.split() if p.startswith("gene:")), "")
            if not gid:
                # ITAG-style: Solyc01g005000.3.1 -> gene Solyc01g005000
                gid = first.split(".")[0]
            name = f"{first}|{gid}"
            buf = []
        elif name:
            buf.append(line.strip().replace("*", ""))
    if name:
        seqs[name] = "".join(buf)
    return seqs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--preset", choices=sorted(PRESETS))
    ap.add_argument("--queries", help="acc=label,acc=label")
    ap.add_argument("--proteome", help="URL of a gzipped protein FASTA")
    ap.add_argument("--label", default="target", help="target species name")
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--all", action="store_true",
                    help="run every job in data/ortholog_jobs.json")
    ap.add_argument("--job", help="run one job from data/ortholog_jobs.json by name")
    a = ap.parse_args()

    if a.all or a.job:
        book = json.loads((DATA / "ortholog_jobs.json").read_text(encoding="utf-8"))["jobs"]
        names = [a.job] if a.job else list(book)
        rc = 0
        for name in names:
            job = book.get(name)
            if not job:
                print(f"no such job: {name}", file=sys.stderr)
                rc = 1
                continue
            print(f"\n########## {name}  ({', '.join(job.get('wanted', []))})")
            rc |= run_job({**job, "entry": name}, a.top)
        return rc

    if a.preset:
        cfg = PRESETS[a.preset]
    elif a.queries and a.proteome:
        cfg = {"label": a.label, "proteome": a.proteome,
               "queries": dict(kv.split("=", 1) for kv in a.queries.split(",")),
               "of_interest": (), "entry": ""}
    else:
        ap.print_help()
        return 1

    return run_job(cfg, a.top)


def run_job(cfg: dict, top: int) -> int:
    try:
        from Bio.Align import PairwiseAligner, substitution_matrices
    except ImportError:
        print("needs biopython: pip install biopython", file=sys.stderr)
        return 1

    print(f"-- queries: {len(cfg['queries'])}")
    qs = query_sequences(cfg["queries"])
    if not qs:
        return 1

    print(f"-- {cfg['label']} proteome: downloading")
    targets = proteome(cfg["proteome"])
    print(f"   {len(targets)} proteins")

    al = PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score, al.mode = -11, -1, "local"

    t0 = time.time()
    best: dict[str, dict[str, float]] = {}
    for i, (tname, tseq) in enumerate(targets.items(), 1):
        if not tseq or len(tseq) > 2000:
            continue
        for label, (_acc, qseq) in qs.items():
            try:
                score = al.score(qseq, tseq)
            except Exception as exc:
                print(f"  ! alignment failed for {label}: {type(exc).__name__}: {exc}",
                      file=sys.stderr)
                return 1
            best.setdefault(tname, {})[label] = float(score)
        if i % 4000 == 0:
            print(f"   {i}/{len(targets)} … {time.time() - t0:.0f}s")

    if not best:
        print("  ! nothing scored — no queries or no target proteins", file=sys.stderr)
        return 1
    rows = []
    for tname, scores in best.items():
        top_label = max(scores, key=scores.get)
        rows.append({
            "protein": tname.split("|")[0],
            "gene": tname.split("|")[1] if "|" in tname else "",
            "closest_query": top_label,
            "scores": {k: round(v, 1) for k, v in sorted(scores.items(), key=lambda x: -x[1])},
            "top_score": round(scores[top_label], 1),
        })
    rows.sort(key=lambda r: -r["top_score"])

    interest = cfg.get("of_interest") or ()
    print(f"\n   {'gene':<24}{'closest query':<46}{'score':>7}")
    shown = 0
    for r in rows:
        if shown >= top:
            break
        mark = "*" if any(r["closest_query"].startswith(x) for x in interest) else " "
        print(f" {mark} {r['gene']:<24}{r['closest_query'][:44]:<46}{r['top_score']:>7}")
        shown += 1
    print("\n   * = closest to the gene of interest rather than to a paralogue")

    payload = {}
    if OUT.exists():
        try:
            payload = json.loads(OUT.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            payload = {}
    payload.setdefault("_about", (
        "Candidate counterparts found by protein similarity with scripts/find_ortholog.py. "
        "Queries include the gene's close paralogues on purpose, so a hit nearer a paralogue is "
        "visible as such. These are ranked candidates, NOT established orthologs: that needs a "
        "phylogenetic tree over the conserved domains and, ideally, expression in the right tissue."
    ))
    payload[cfg.get("entry") or cfg["label"]] = {
        "target_species": cfg["label"],
        "proteome": cfg["proteome"],
        "queries": {label: QUERY_TITLES.get(label, acc) for label, (acc, _s) in qs.items()},
        "checked_on": time.strftime("%Y-%m-%d"),
        "candidates": rows[:top],
    }
    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"\nwrote {OUT.relative_to(ROOT)} in {time.time() - t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
