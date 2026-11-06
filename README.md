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

1. **Mine** functional-evidence candidate genes per trait from PubMed (`scripts/mine_candidates.py`, run via Actions). The query *requires* experimental language (overexpression, knockout, CRISPR, mutant…) so associations-only hits are excluded.
2. **Curate** candidates into the schema (`data/atlas.json`) — a human (breeder) decides what is real and fills portability, linkage, side-effects.
3. **IP layer** (planned): fill the patent signal automatically from open patent data (PatentsView / Lens).
4. **Publish** the site from `data/atlas.json`.

## Status

Proof-of-concept seed: a few curated, functionally-validated entries for **salt tolerance** and **seedlessness** (`data/atlas.json`). Not a complete atlas.

## Honesty

Hypotheses and research signals, not verdicts. AI and open data point to candidates; the breeder's judgement decides.
