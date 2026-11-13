# Transfer Atlas

*An open map of genetic solutions that already exist — so breeders stop solving the same trait from scratch.*

Part of [AgriHub-AI](https://agrihubai.wordpress.com) · Almería, 2026.

## The idea

For a given **trait** and **target crop**, the Atlas shows where a genetic solution is *already documented* — in another crop or a wild relative — and tags how reachable it is in the gene-editing / NGT era. The EU adopted its New Genomic Techniques (NGT) framework in 2026 (applying from ~2028), so breeders now need a target list; this is that list.

Your sentence, if you are a breeder: *"I want trait X in crop Y."* The Atlas answers: *who already solved it, with which gene, how portable it is to your crop, what baggage it carries, and whether the space is already patented.*

## What each entry carries

| Field | Meaning | Kind |
|---|---|---|
| gene / source / does | the solution and where it is documented | **sourced fact** |
| evidence (+ link) | functional proof (overexpression, knockout, field…) | **sourced fact** |
| target note | status in your crop (present / variant / missing) | indicative |
| **portability** | native edit (NGT cat-1) · cisgenic · transgenic | indicative |
| **linkage** | linkage drag if introgressed — and how native editing avoids it | indicative |
| **side-effects** | known pleiotropy / trade-offs | indicative |
| **IP** | patent/application activity on the gene or use | research signal, **not legal advice** |

Gene, source and evidence are sourced. Portability, linkage, side-effects and IP are **indicative research signals** to guide a breeder's own check — not lab results or legal classifications. Every candidate still needs allele checking and validation in the target crop.

## How it is built

1. **Curate** entries into `data/atlas_v2.json` (the draft, left untouched for provenance) and
   record every curated decision in `data/corrections.json` — who decided, why, and from what.
2. **Verify** with `scripts/verify_genes.py`, which asks NCBI Gene, PubMed, Sol Genomics Network
   and the Ensembl Plants annotation what they actually hold, and writes only what they return.
   Anything it cannot confirm stays empty and the entry says so.
3. **Resolve locus tags** with `scripts/resolve_loci_ensembl.py`, which streams Ensembl Plants
   GFF3 files and writes `data/locus_resolutions.json`, so locus checks need no network and give
   the same answer every run.
4. **Publish** the site from `data/atlas_v3.json`.

Run it yourself:

```bash
python3 scripts/verify_genes.py --migrate --verify --accessions
python3 scripts/resolve_loci_ensembl.py
```

## Status

31 entries across 10 vegetable crops, as of 13 November 2026:

| | entries |
|---|---|
| **verified** — every named gene resolves, and a reference is confirmed | 18 |
| **partial** — some genes or some references resolve | 12 |
| **unverified** — neither resolved | 1 |

Evidence level is **derived from what was verified**, never from what a draft claimed: 2 at
Level 1, 16 at Level 2, 7 at Level 3, 4 at Level 4, 2 at Level 5. Level 1 is never derived
automatically — it needs a person to name a commercial cultivar, and the product has to embody
*this* entry's genetic solution rather than the same trait by another route.

Locus identifiers: 24 of 28 confirmed. The four that are not need the eggplant and watermelon
genome databases, which are unreachable from the build environment.

`docs/verification-gaps.md` lists what is still open and what would close it.

## What verification means here

A reference counts only when PubMed matches the cited author, year and journal, or when a
supplied PMID resolves. A PMID found by a topic search alone does **not** count: it may be a
different paper about the same gene, and it is labelled as unconfirmed in the data.

Gene identity comes from an NCBI Gene record, a GenBank record carrying the locus tag, or a
Sol Genomics / Ensembl Plants annotation. For several crops — eggplant, onion, watermelon — no
NCBI gene build exists, so an accession is the anchor.

## Honesty

Gene identity and references are machine-verified and re-runnable. **Portability, linkage drag,
pleiotropy and IP remain indicative research signals curated by a breeder** — not lab results
and not legal classifications. Every candidate still needs allele checking and validation in the
target crop.

Where something is wrong, the file says so rather than dropping it quietly: `corrections.json`
records superseded genes, rejected commercial claims, identifiers that turned out to belong to
other organisms, and open questions about identifiers that do not match their entry.
