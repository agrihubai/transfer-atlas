# What is still open, and what would close it

Generated from `data/atlas_v3.json` on 16 November 2026 by `scripts/write_gaps.py`. Do not edit by hand — rerun the script.

State: **22 verified, 10 partial, 0 unverified** of 32 entries.

## What counts, so nobody wastes effort

A reference counts only when PubMed matches the cited **author, year and journal**, or when a supplied PMID resolves. A PMID found by topic search alone does not count: it may be a different paper about the same gene.

**The most useful thing anyone can send** is a gene symbol, its crop species, and a DOI or PMID for a paper that functionally tests that gene in that crop. The pipeline resolves the rest. Reconstructed citations do not help — an author and year with no findable journal is what the original draft had, and it is what produced most of these gaps.

**Not genes, so do not look for them:** `uORF` (a 5' UTR feature), `ROS scavengers`, `ABA catabolism`, `prickle loci`. These were descriptive phrases in the draft's gene field that the parser read as symbols.

## Already searched — please do not repeat these

Every unresolved symbol below has been looked for, by script, in all of these:

- **NCBI Gene**, as written and with the species prefix stripped (`SlGAD3` and `GAD3`), constrained to the stated organism, with hits ranked by whether the symbol is the official name or only an alias.
- **NCBI nuccore and protein**, by locus tag, constrained by organism.
- **Sol Genomics Network** feature search, for Solanaceae symbols and locus tags.
- **Ensembl Plants annotation files** (ftp.ebi.ac.uk), streamed per species and searched both by locus tag and by gene name — cucumber, melon, watermelon, both Brassicas, pepper and carrot. Cucumber tags also tried in Ensembl's underscore form (`Csa_5G623470`).
- **Eggplant V4 annotation** (hosted by Sol Genomics), including the V3-to-V4 conversion table.

These annotations mostly carry no gene names, which is why a name search found almost nothing: of 21 unresolved symbols it resolved none. So the useful answer is **an accession or gene ID read off a record or a paper's data-availability section** — not another database name, and not a guess.

**What is unreachable from the build environment**, so do not suggest it: Ensembl's REST service for plants (rest.ensemblgenomes.org serves a certificate for another host), plants.ensembl.org (redirects to an unreachable host), cucurbitgenomics.org (no working HTTPS), and the watermelon 97103 v2 build (nowhere reachable).

**Three rules, from what went wrong last time.** Identifiers supplied earlier included a fungal ITS sequence given as a watermelon gene, a bacterial gene and a tomato protein given as eggplant polyphenol oxidases, an NCBI id for a protein kinase given as a lettuce GGP2, and two PMIDs that were a cellulose-microfibril paper and a cardiovascular review. The PMIDs supplied were otherwise reliable. So:

1. Never give an identifier not read off a record. Write **unknown** instead — an honest gap is useful, a wrong accession costs more than it saves.
2. Always name the **journal** with author and year. The matcher requires it, because author and year alone cannot separate two papers.
3. Say what a paper does **not** cover. That is what saved several entries from false claims.

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

## Genes resolve, no reference is confirmed (4)

Identity is settled; the functional evidence is not. A DOI or PMID for a knockout, overexpression, allele-swap or field study closes these.

- **CUC-SEX-01** — Cucumber, Gynoecy & parthenocarpy (Female flower dominance & seedless fruit) (`CsACS2, CsACS11, CsWIP1, CsARF8`)
  - cited: "Hu et al., 2017, Plant J." — not matched
- **CAR-COLOR-01** — Carrot, Taproot carotenoid & beta-carotene accumulation (`DcRPGE1, DcAPRR2, DcPSY1, DcPSY2, DcLCYB1`)
  - cited: "2024 Functional Carrot Root Color Study; Ellison et al., 2018" — not matched
- **GAP-RECOV-01** — Tomato, Post-stress growth & yield recovery kinetics (`CYP707A1`)
  - cited: "2025–2026 Plant Resilience Frameworks" — not matched
- **WAT-SEX-01** — Watermelon, Gynoecy / sex expression (`ClWIP1`)
  - cited: "" — not matched

## Nothing resolves: neither gene nor reference (1)

These should be filled with a real source or retired from the atlas.

- **CUC-MILDEW-01** — Cucumber, Powdery mildew resistance (`CsMLO1, CsMLO8, CsMLO11`)
  - cited: "Berg et al., 2015, Front. Plant Sci."

## Some genes of an entry resolve, others do not (7)

- **TOM-DWARF-01** — Tomato (3/5 resolved): `SlER`, `SlYTH2`
- **TOM-SIZE-01** — Tomato (1/3 resolved): `SlKIX8`, `SlKIX9`
- **CUC-SEX-01** — Cucumber (2/4 resolved): `CsWIP1`, `CsARF8`
- **MEL-SHELF-01** — Melon (2/3 resolved): `CmNAC-NOR`
- **LET-BOLT-01** — Lettuce (1/4 resolved): `LsCO`, `LsRGL1`, `LsSOC1`
- **BRA-GLUC-01** — Brassica (Broccoli/Cabbage/Cauliflower) (1/4 resolved): `HAG1`, `MYB29`, `MYB34`
- **CAR-COLOR-01** — Carrot (3/5 resolved): `DcRPGE1`, `DcAPRR2`

## Every gene symbol with no identifier, by species (18)

For each: the entry, the symbol as the atlas writes it, and the name that was searched after stripping the species prefix. An accession, a gene ID, or `unknown`.

**Brassica oleracea**

- `HAG1` (BRA-GLUC-01) — searched as `HAG1`
- `MYB29` (BRA-GLUC-01) — searched as `MYB29`
- `MYB34` (BRA-GLUC-01) — searched as `MYB34`

**Cucumis melo**

- `CmNAC-NOR` (MEL-SHELF-01) — searched as `NAC-NOR`

**Cucumis sativus**

- `CsMLO1` (CUC-MILDEW-01) — searched as `MLO1`
- `CsMLO8` (CUC-MILDEW-01) — searched as `MLO8`
- `CsMLO11` (CUC-MILDEW-01) — searched as `MLO11`
- `CsWIP1` (CUC-SEX-01) — searched as `WIP1`
- `CsARF8` (CUC-SEX-01) — searched as `ARF8`

**Daucus carota**

- `DcRPGE1` (CAR-COLOR-01) — searched as `RPGE1`
- `DcAPRR2` (CAR-COLOR-01) — searched as `APRR2`

**Lactuca sativa**

- `LsCO` (LET-BOLT-01) — searched as `CO`
- `LsRGL1` (LET-BOLT-01) — searched as `RGL1`
- `LsSOC1` (LET-BOLT-01) — searched as `SOC1`

**Solanum lycopersicum**

- `SlER` (TOM-DWARF-01) — searched as `ER`
- `SlYTH2` (TOM-DWARF-01) — searched as `YTH2`
- `SlKIX8` (TOM-SIZE-01) — searched as `KIX8`
- `SlKIX9` (TOM-SIZE-01) — searched as `KIX9`

## Locus tags that no reachable database carries (2)

Not evidence that a tag is wrong — the database that holds it is unreachable, or the tag comes from a different genome build than the one available. **The useful answer here is which genome build each tag belongs to**, so it can be converted.

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
