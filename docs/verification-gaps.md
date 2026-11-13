# What is still open, and what would close it

Generated from `data/atlas_v3.json` on 13 November 2026 by `scripts/write_gaps.py`. Do not edit by hand — rerun the script.

State: **18 verified, 12 partial, 1 unverified** of 31 entries.

## What counts, so nobody wastes effort

A reference counts only when PubMed matches the cited **author, year and journal**, or when a supplied PMID resolves. A PMID found by topic search alone does not count: it may be a different paper about the same gene.

**The most useful thing anyone can send** is a gene symbol, its crop species, and a DOI or PMID for a paper that functionally tests that gene in that crop. The pipeline resolves the rest. Reconstructed citations do not help — an author and year with no findable journal is what the original draft had, and it is what produced most of these gaps.

**Not genes, so do not look for them:** `uORF` (a 5' UTR feature), `ROS scavengers`, `ABA catabolism`, `prickle loci`. These were descriptive phrases in the draft's gene field that the parser read as symbols.

## Needs one look, nothing more (2)

The author and journal match the citation, but the title does not carry the trait, so the pipeline will not count it. A person can settle this in seconds.

- **TOM-FIRM-01** — Tomato, Fruit firmness & extended shelf life without ripening block
  - cited: Uluisik et al., 2016, Nat. Biotechnol.; Wang et al., 2019
  - candidate: PMID 27454737 (Nat Biotechnol) Genetic improvement of tomato by targeted control of fruit softening.
  - question: is this the intended paper? yes or no
- **GAP-COMB-01** — Tomato, Combined Heat + Salinity + Drought Tolerance (Mediterranean field conditions)
  - cited: Zandalinas et al., 2021; Rivero et al., 2022
  - candidate: PMID 33496342 (New Phytol) The impact of multifactorial stress combination on plant growth and survival.
  - question: is this the intended paper? yes or no

## A confirmed paper exists, but no gene symbol resolves (2)

The cheapest wins: the evidence is confirmed, only the identifier is missing. The paper's own data-availability section usually has it.

- **EGG-BROWN-01** — Eggplant (Solanum melongena), Non-browning cut fruit flesh
  - unresolved: `SmelPPO4`, `SmelPPO5`, `SmelPPO6`
  - confirmed paper: PMID 33343607 — Simultaneous CRISPR/Cas9 Editing of Three PPO Genes Reduces Fruit Flesh Browning in Sola
  - needed: an accession or gene ID for those symbols, read off a record
- **EGG-PRICKLE-01** — Eggplant (Solanum melongena), Prickle-free calyx & stem
  - unresolved: `SmLOG1`
  - confirmed paper: PMID 38974191 — PE (Prickly Eggplant) encoding a cytokinin-activating enzyme responsible for the formati
  - needed: an accession or gene ID for those symbols, read off a record

## Genes resolve, no reference is confirmed (4)

Identity is settled; the functional evidence is not. A DOI or PMID for a knockout, overexpression, allele-swap or field study closes these.

- **CUC-SEX-01** — Cucumber, Gynoecy & parthenocarpy (Female flower dominance & seedless fruit) (`CsACS2, CsACS11, CsWIP1, CsARF8`)
  - cited: "Hu et al., 2017, Plant J." — not matched
- **CAR-COLOR-01** — Carrot, Taproot carotenoid & beta-carotene accumulation (`DcRPGE1, DcAPRR2, DcPSY1, DcPSY2, DcLCYB1`)
  - cited: "2024 Functional Carrot Root Color Study; Ellison et al., 2018" — not matched
- **GAP-COMB-01** — Tomato, Combined Heat + Salinity + Drought Tolerance (Mediterranean field conditions) (`SOS1, SOS2, NHX1, DREB2, HsfA2`)
  - cited: "Zandalinas et al., 2021; Rivero et al., 2022" — not matched
- **GAP-RECOV-01** — Tomato, Post-stress growth & yield recovery kinetics (`CYP707A1`)
  - cited: "2025–2026 Plant Resilience Frameworks" — not matched

## Nothing resolves: neither gene nor reference (2)

These should be filled with a real source or retired from the atlas.

- **CUC-MILDEW-01** — Cucumber, Powdery mildew resistance (`CsMLO1, CsMLO8, CsMLO11`)
  - cited: "Berg et al., 2015, Front. Plant Sci."
- **WAT-SEED-01** — Watermelon, Diploid seedless / stenospermocarpic fruit development (`ClAGL11 (plus ClWIP1 for gynoecy)`)
  - cited: "Wang et al., 2021; Cucurbit AGL11 studies"

## Some genes of an entry resolve, others do not (7)

- **TOM-DWARF-01** — Tomato (3/5 resolved): `SlER`, `SlYTH2`
- **TOM-SIZE-01** — Tomato (1/3 resolved): `SlKIX8`, `SlKIX9`
- **CUC-SEX-01** — Cucumber (2/4 resolved): `CsWIP1`, `CsARF8`
- **MEL-SHELF-01** — Melon (2/3 resolved): `CmNAC-NOR`
- **LET-BOLT-01** — Lettuce (1/4 resolved): `LsCO`, `LsRGL1`, `LsSOC1`
- **BRA-GLUC-01** — Brassica (Broccoli/Cabbage/Cauliflower) (1/4 resolved): `HAG1`, `MYB29`, `MYB34`
- **CAR-COLOR-01** — Carrot (3/5 resolved): `DcRPGE1`, `DcAPRR2`

## Locus tags that no reachable database carries (4)

Not evidence that a tag is wrong — the database that holds it is unreachable, or the tag comes from a different genome build than the one available.

- `Cla97C07G144940` — WAT-SEED-01, Citrullus lanatus
- `SMEL4.1_08g017350` — EGG-BROWN-01, Solanum melongena
- `SMEL4.1_06g...` — EGG-PRICKLE-01, Solanum melongena
- `BraA03g012133E` — BRA-CLUB-01, Brassica rapa

## Commercial status, which only a person can settle (6)

Level 1 requires a named cultivar or product whose trait comes from **this entry's genetic solution**, not the same trait by another route.

- **TOM-TYLCV-01** — Tomato, TYLCV resistance (ty-5 locus)
  - NOT set to Level 1. The cultivar names offered — 'Camaro', 'Grand Marshall', 'Jolene' — were not confirmed, and the claim mixed introgressed natural ty-5 alleles with targeted modifications, which are different routes. ty-5 resistance is wi
- **PEP-PUNG-01** — Pepper, Pungency control (Capsaicinoid presence vs. non-pungent)
  - NOT set to Level 1 by the pipeline, and worth a human decision. Every commercial sweet bell pepper is pun1 loss-of-function, so the allele is commercially ubiquitous — but it arrived through domestication centuries ago, not through this sol
- **PEP-VIRUS-01** — Pepper, Potyvirus resistance (PVY, TEV, PepMoV — pvr1 / pvr2)
  - NOT set to Level 1: no cultivar was named, only 'commercial hybrid lines bred with pvr1/pvr2' and a list of seed companies. A category is not a product.
- **CUC-BITTER-01** — Cucumber, Non-bitter fruit while retaining bitter-leaf defense
  - NOT set to Level 1. 'Marketmore' and 'Spacemaster' are long-standing cucumber cultivars, and the claim that they were 'engineered or selected' for tissue-specific bitterness regulation conflates ordinary non-bitter fruit with the Bt/Bl regu
- **BRA-CLUB-01** — Brassica rapa, Clubroot resistance (Plasmodiophora brassicae)
  - NOT set to Level 1. The literature does document that CR loci have been introduced into commercial F1 cultivars of Chinese cabbage, and that CRb/CRbKato is the same allele as CRa (Breeding Science 72:21040; Hatakeyama et al. 2017), and Syng
- **ONI-TEAR-01** — Onion, Tearless / low-lachrymatory onion & reduced pungency
  - NOT set to Level 1, and this is the clearest case of why the rule matters. Sunions is a real marketed tearless onion from BASF/Nunhems, but BASF states it is the result of more than three decades of CONVENTIONAL breeding, with pungency trac
