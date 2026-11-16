# What is still open, and what would close it

Generated from `data/atlas_v3.json` on 16 November 2026 by `scripts/write_gaps.py`. Do not edit by hand — rerun the script.

State: **0 verified, 0 partial, 32 unverified** of 32 entries.

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


## Nothing resolves: neither gene nor reference (32)

These should be filled with a real source or retired from the atlas.

- **TOM-GABA-01** — Tomato, High GABA accumulation (`SlGAD3, SlGAD2`)
  - cited: "Nonaka et al., 2017 (Sci. Rep.); 2025 Elite Genotype Evaluation"
- **TOM-VITD-01** — Tomato, Provitamin D3 (7-DHC) biofortification (`Sl7-DR2`)
  - cited: "Li et al., 2022, Nature Plants"
- **TOM-TOBRFV-01** — Tomato, ToBRFV resistance (Tobamovirus multiplication) (`SlTOM1a, SlTOM3`)
  - cited: "Ishikawa et al., 2022; Kravchik et al., 2022"
- **TOM-TYLCV-01** — Tomato, TYLCV resistance (ty-5 locus) (`SlPelo`)
  - cited: "Lapidot et al., 2015; Pramanik et al., 2021"
- **TOM-MILDEW-01** — Tomato, Powdery mildew resistance (Oidium neolycopersici) (`SlMlo1`)
  - cited: "Nekrasov et al., 2017, Sci. Rep."
- **TOM-BROAD-01** — Tomato, Multi-pathogen resistance (Bacterial speck, Phytophthora, Powdery mildew) (`SlDMR6-1`)
  - cited: "Thomazella et al., 2021, PNAS"
- **TOM-DWARF-01** — Tomato, Compact / dwarf architecture for vertical farming & early harvest (`SlER, SP5G, SP, SlYABBY2b, SlYTH2`)
  - cited: "Kwon et al., 2020, Nat. Biotechnol."
- **TOM-PARTH-01** — Tomato, Parthenocarpy (seedless fruit set — distinct from heat tolerance) (`SlIAA9, SlAGL6, SlARF8`)
  - cited: "Ueta et al., 2017; Klap et al., 2017"
- **TOM-SIZE-01** — Tomato, Fruit size & pericarp thickness (`FW2.2, SlKIX8, SlKIX9`)
  - cited: "Swinnen et al., 2022, Plant Physiol."
- **TOM-FIRM-01** — Tomato, Fruit firmness & extended shelf life without ripening block (`SlPL, SlPG2a, SlCEL2`)
  - cited: "Uluisik et al., 2016, Nat. Biotechnol.; Wang et al., 2019"
- **PEP-PUNG-01** — Pepper, Pungency control (Capsaicinoid presence vs. non-pungent) (`Pun1 (AT3)`)
  - cited: "Stewart et al., 2005, Plant J."
- **PEP-PAMT-01** — Pepper, Low-pungency capsinoid accumulation (`pAMT`)
  - cited: "Lang et al., 2009, Plant J."
- **PEP-VIRUS-01** — Pepper, Potyvirus resistance (PVY, TEV, PepMoV — pvr1 / pvr2) (`eIF4E1, eIFiso4E`)
  - cited: "Ruffel et al., 2002; Yoon et al., 2020"
- **PEP-MILDEW-01** — Pepper, Powdery mildew resistance (Leveillula taurica) (`CaMlo2, CaMlo1`)
  - cited: "Zheng et al., 2013; Kim & Hwang, 2012"
- **CUC-BITTER-01** — Cucumber, Non-bitter fruit while retaining bitter-leaf defense (`Bt, Bl, Bi`)
  - cited: "Shang et al., 2014, Science"
- **CUC-MILDEW-01** — Cucumber, Powdery mildew resistance (`CsMLO1, CsMLO8, CsMLO11`)
  - cited: "Berg et al., 2015, Front. Plant Sci."
- **CUC-SEX-01** — Cucumber, Gynoecy & parthenocarpy (Female flower dominance & seedless fruit) (`CsACS2, CsACS11, CsWIP1, CsARF8`)
  - cited: "Hu et al., 2017, Plant J."
- **MEL-SHELF-01** — Melon, Extended shelf life & climacteric ripening control (`CmACO1, CmACS7, CmNAC-NOR`)
  - cited: "Nonaka et al., 2023; Ayub et al., 1996"
- **WAT-SEED-01** — Watermelon, Seedlessness (stenospermocarpy) — documented in grape, untested in watermelon (`AGL11`)
  - cited: "Wang et al., 2021; Cucurbit AGL11 studies"
- **WAT-BRIX-01** — Watermelon, Sugar accumulation (Brix) (`ClTST2`)
  - cited: "Ren et al., 2018, Plant Cell; Guo et al., 2019"
- **EGG-BROWN-01** — Eggplant, Non-browning cut fruit flesh (`SmelPPO4, SmelPPO5, SmelPPO6`)
  - cited: "Maioli et al., 2020, Front. Plant Sci."
- **EGG-PRICKLE-01** — Eggplant, Prickle-free calyx & stem (`SmLOG1`)
  - cited: "Miyatake et al., 2020; Solanaceae prickle genetics"
- **LET-BOLT-01** — Lettuce, Delayed bolting (premature flowering under summer heat) (`LsFT, LsCO, LsRGL1, LsSOC1`)
  - cited: "Chen et al., 2018, Plant Biotechnol. J."
- **LET-VITC-01** — Lettuce, Elevated Vitamin C (Ascorbate) via uORF editing (`LsGGP2`)
  - cited: "Zhang et al., 2018, Nat. Biotechnol."
- **BRA-CLUB-01** — Brassica rapa, Clubroot resistance (Plasmodiophora brassicae) (`CRa`)
  - cited: "Ueno et al., 2012; Hatakeyama et al., 2013"
- **BRA-GLUC-01** — Brassica (Broccoli/Cabbage/Cauliflower), Glucoraphanin / Sulforaphane biofortification (`MYB28 (HAG1), MYB29, MYB34`)
  - cited: "Traka et al., 2013, New Phytol."
- **CAR-COLOR-01** — Carrot, Taproot carotenoid & beta-carotene accumulation (`DcRPGE1, DcAPRR2, DcPSY1, DcPSY2, DcLCYB1`)
  - cited: "2024 Functional Carrot Root Color Study; Ellison et al., 2018"
- **ONI-TEAR-01** — Onion, Tearless / low-lachrymatory onion & reduced pungency (`LFS`)
  - cited: "Eady et al., 2008, Plant Physiol."
- **GAP-COMB-01** — Tomato, Combined Heat + Salinity + Drought Tolerance (Mediterranean field conditions) (`SOS1, SOS2, NHX1, DREB2, HsfA2`)
  - cited: "Zandalinas et al., 2021; Rivero et al., 2022"
- **GAP-REPRO-01** — Tomato, Pollen & ovule thermotolerance + stigma receptivity under heat waves (`HsfA2, HSP70, APX1`)
  - cited: "Rieu et al., 2017; Fragkostefanakis et al., 2016"
- **GAP-RECOV-01** — Tomato, Post-stress growth & yield recovery kinetics (`CYP707A1`)
  - cited: "2025–2026 Plant Resilience Frameworks"
- **WAT-SEX-01** — Watermelon, Gynoecy / sex expression (`ClWIP1`)
  - cited: ""

## Every gene symbol with no identifier, by species (74)

For each: the entry, the symbol as the atlas writes it, and the name that was searched after stripping the species prefix. An accession, a gene ID, or `unknown`.

**Allium cepa**

- `LFS` (ONI-TEAR-01) — searched as `LFS`

**Brassica oleracea**

- `MYB28` (BRA-GLUC-01) — searched as `MYB28`
- `HAG1` (BRA-GLUC-01) — searched as `HAG1`
- `MYB29` (BRA-GLUC-01) — searched as `MYB29`
- `MYB34` (BRA-GLUC-01) — searched as `MYB34`

**Brassica rapa**

- `CRa` (BRA-CLUB-01) — searched as `CRa`

**Capsicum annuum**

- `Pun1` (PEP-PUNG-01) — searched as `Pun1`
- `AT3` (PEP-PUNG-01) — searched as `AT3`
- `pAMT` (PEP-PAMT-01) — searched as `pAMT`
- `eIF4E1` (PEP-VIRUS-01) — searched as `eIF4E1`
- `eIFiso4E` (PEP-VIRUS-01) — searched as `eIFiso4E`
- `CaMlo2` (PEP-MILDEW-01) — searched as `Mlo2`
- `CaMlo1` (PEP-MILDEW-01) — searched as `Mlo1`

**Citrullus lanatus**

- `AGL11` (WAT-SEED-01) — searched as `AGL11`
- `ClTST2` (WAT-BRIX-01) — searched as `TST2`
- `ClWIP1` (WAT-SEX-01) — searched as `WIP1`

**Cucumis melo**

- `CmACO1` (MEL-SHELF-01) — searched as `ACO1`
- `CmACS7` (MEL-SHELF-01) — searched as `ACS7`
- `CmNAC-NOR` (MEL-SHELF-01) — searched as `NAC-NOR`

**Cucumis sativus**

- `Bt` (CUC-BITTER-01) — searched as `Bt`
- `Bl` (CUC-BITTER-01) — searched as `Bl`
- `Bi` (CUC-BITTER-01) — searched as `Bi`
- `CsMLO1` (CUC-MILDEW-01) — searched as `MLO1`
- `CsMLO8` (CUC-MILDEW-01) — searched as `MLO8`
- `CsMLO11` (CUC-MILDEW-01) — searched as `MLO11`
- `CsACS2` (CUC-SEX-01) — searched as `ACS2`
- `CsACS11` (CUC-SEX-01) — searched as `ACS11`
- `CsWIP1` (CUC-SEX-01) — searched as `WIP1`
- `CsARF8` (CUC-SEX-01) — searched as `ARF8`

**Daucus carota**

- `DcRPGE1` (CAR-COLOR-01) — searched as `RPGE1`
- `DcAPRR2` (CAR-COLOR-01) — searched as `APRR2`
- `DcPSY1` (CAR-COLOR-01) — searched as `PSY1`
- `DcPSY2` (CAR-COLOR-01) — searched as `PSY2`
- `DcLCYB1` (CAR-COLOR-01) — searched as `LCYB1`

**Lactuca sativa**

- `LsFT` (LET-BOLT-01) — searched as `FT`
- `LsCO` (LET-BOLT-01) — searched as `CO`
- `LsRGL1` (LET-BOLT-01) — searched as `RGL1`
- `LsSOC1` (LET-BOLT-01) — searched as `SOC1`
- `LsGGP2` (LET-VITC-01) — searched as `GGP2`

**Solanum lycopersicum**

- `SlGAD3` (TOM-GABA-01) — searched as `GAD3`
- `SlGAD2` (TOM-GABA-01) — searched as `GAD2`
- `Sl7-DR2` (TOM-VITD-01) — searched as `7-DR2`
- `SlTOM1a` (TOM-TOBRFV-01) — searched as `TOM1a`
- `SlTOM3` (TOM-TOBRFV-01) — searched as `TOM3`
- `SlPelo` (TOM-TYLCV-01) — searched as `Pelo`
- `SlMlo1` (TOM-MILDEW-01) — searched as `Mlo1`
- `SlDMR6-1` (TOM-BROAD-01) — searched as `DMR6-1`
- `SlER` (TOM-DWARF-01) — searched as `ER`
- `SP5G` (TOM-DWARF-01) — searched as `SP5G`
- `SP` (TOM-DWARF-01) — searched as `SP`
- `SlYABBY2b` (TOM-DWARF-01) — searched as `YABBY2b`
- `SlYTH2` (TOM-DWARF-01) — searched as `YTH2`
- `SlIAA9` (TOM-PARTH-01) — searched as `IAA9`
- `SlAGL6` (TOM-PARTH-01) — searched as `AGL6`
- `SlARF8` (TOM-PARTH-01) — searched as `ARF8`
- `FW2.2` (TOM-SIZE-01) — searched as `FW2.2`
- `SlKIX8` (TOM-SIZE-01) — searched as `KIX8`
- `SlKIX9` (TOM-SIZE-01) — searched as `KIX9`
- `SlPL` (TOM-FIRM-01) — searched as `PL`
- `SlPG2a` (TOM-FIRM-01) — searched as `PG2a`
- `SlCEL2` (TOM-FIRM-01) — searched as `CEL2`
- `SOS1` (GAP-COMB-01) — searched as `SOS1`
- `SOS2` (GAP-COMB-01) — searched as `SOS2`
- `NHX1` (GAP-COMB-01) — searched as `NHX1`
- `DREB2` (GAP-COMB-01) — searched as `DREB2`
- `HsfA2` (GAP-COMB-01) — searched as `HsfA2`
- `HsfA2` (GAP-REPRO-01) — searched as `HsfA2`
- `HSP70` (GAP-REPRO-01) — searched as `HSP70`
- `APX1` (GAP-REPRO-01) — searched as `APX1`
- `CYP707A1` (GAP-RECOV-01) — searched as `CYP707A1`

**Solanum melongena**

- `SmelPPO4` (EGG-BROWN-01) — searched as `PPO4`
- `SmelPPO5` (EGG-BROWN-01) — searched as `PPO5`
- `SmelPPO6` (EGG-BROWN-01) — searched as `PPO6`
- `SmLOG1` (EGG-PRICKLE-01) — searched as `LOG1`

## Commercial status, which only a person can settle (6)

Level 1 requires a named cultivar or product whose trait comes from **this entry's genetic solution**, not the same trait by another route.

- **TOM-TYLCV-01** — Tomato, TYLCV resistance (ty-5 locus)
  - NOT set to Level 1. The cultivar names offered — 'Camaro', 'Grand Marshall', 'Jolene' — were not confirmed, and the claim mixed introgressed natural ty-5 alleles with targeted modifications, which are different routes. ty-5 resistance is wi
- **PEP-PUNG-01** — Pepper, Pungency control (Capsaicinoid presence vs. non-pungent)
  - NOT Level 1 under the rule settled on 2026-11-16 (see _about.level1_rule). Every commercial sweet bell pepper is pun1 loss-of-function, so the allele is commercially ubiquitous — but it arrived through domestication centuries ago, not by an
- **PEP-VIRUS-01** — Pepper, Potyvirus resistance (PVY, TEV, PepMoV — pvr1 / pvr2)
  - NOT set to Level 1: no cultivar was named, only 'commercial hybrid lines bred with pvr1/pvr2' and a list of seed companies. A category is not a product.
- **CUC-BITTER-01** — Cucumber, Non-bitter fruit while retaining bitter-leaf defense
  - NOT set to Level 1. 'Marketmore' and 'Spacemaster' are long-standing cucumber cultivars, and the claim that they were 'engineered or selected' for tissue-specific bitterness regulation conflates ordinary non-bitter fruit with the Bt/Bl regu
- **BRA-CLUB-01** — Brassica rapa, Clubroot resistance (Plasmodiophora brassicae)
  - NOT set to Level 1. The literature does document that CR loci have been introduced into commercial F1 cultivars of Chinese cabbage, and that CRb/CRbKato is the same allele as CRa (Breeding Science 72:21040; Hatakeyama et al. 2017), and Syng
- **ONI-TEAR-01** — Onion, Tearless / low-lachrymatory onion & reduced pungency
  - NOT set to Level 1, and this is the clearest case of why the rule matters. Sunions is a real marketed tearless onion from BASF/Nunhems, but BASF states it is the result of more than three decades of CONVENTIONAL breeding, with pungency trac
