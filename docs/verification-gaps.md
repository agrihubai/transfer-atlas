# Still unverified — what a second opinion could actually fix

State on 2026-11-13 after running `scripts/verify_genes.py` against NCBI Gene and PubMed, applying the maintainer's curated corrections, and deriving evidence levels from what was verified.

## Read this first

**A reference counts only if PubMed matched the cited author, year and gene, or if a supplied PMID resolved.** A PMID found by topic search alone does not count — it may be a different paper about the same gene. Those are listed below as unconfirmed.

**What would help, concretely:** gene symbol + crop species + a DOI or PMID of a paper that functionally tests that gene in that crop. The pipeline resolves everything else. Please do not send reconstructed citations — an author and year with no findable DOI is what the draft already had, and it is what produced most of the gaps below.

**Do not chase these, they are not genes:** `uORF` (a feature of LsGGP2), `ROS scavengers`, `ABA catabolism`, `for gynoecy`, `prickle loci`. Descriptive phrases in the draft's gene field were parsed as symbols. Likewise these are one gene under two names: `PLETHORA`/`PLT`, `HAG1`/`MYB28`, `AT3`/`Pun1`, and `Bt`/`CsBt`, `Bl`/`CsBl`, `Bi`/`CsBi` (Cs is the species prefix).

**Locus IDs are a separate problem, already diagnosed — skip them.** Genome-project tags (Solyc, Csa, Cla97, MELO3C, SMEL4.1, BraA) are not indexed in NCBI, and the plant databases that hold them are unreachable from the build environment. 25 of 28 are unresolved for that reason alone, not because they look wrong.

## Priority 1 — a confirmed paper exists, but no gene symbol resolves (11 entries)

These are the cheapest wins: the functional evidence is already confirmed in PubMed. What is missing is the gene symbol or identifier that the paper itself uses, in the stated species. Reading the paper's methods usually settles it.

- **TOM-VITD-01** — Tomato, Provitamin D3 (7-DHC) biofortification
  - symbols that do not resolve: `Sl7-DR2`
  - confirmed paper: PMID 35606499 — Biofortified tomatoes provide a new route to vitamin D sufficiency.
  - needed: the gene symbol used in that paper for Solanum lycopersicum, or an NCBI/Ensembl gene ID
- **TOM-TOBRFV-01** — Tomato, ToBRFV resistance (Tobamovirus multiplication)
  - symbols that do not resolve: `SlTOM1a`, `SlTOM3`
  - confirmed paper: PMID 35706371 — Knockout of SlTOM1 and SlTOM3 results in differential resistance to tobamovirus in tomato.
  - needed: the gene symbol used in that paper for Solanum lycopersicum, or an NCBI/Ensembl gene ID
- **TOM-TYLCV-01** — Tomato, TYLCV resistance (ty-5 locus)
  - symbols that do not resolve: `SlPelo`
  - confirmed paper: PMID 26448569 — A Novel Route Controlling Begomovirus Resistance by the Messenger RNA Surveillance Factor Pelot
  - needed: the gene symbol used in that paper for Solanum lycopersicum, or an NCBI/Ensembl gene ID
- **PEP-MILDEW-01** — Pepper, Powdery mildew resistance (Leveillula taurica)
  - symbols that do not resolve: `CaMLO2 (CaMLO`
  - confirmed paper: PMID 23923019 — Loss of function in Mlo orthologs reduces susceptibility of pepper and tomato to powdery mildew
  - needed: the gene symbol used in that paper for Capsicum annuum, or an NCBI/Ensembl gene ID
- **CUC-BITTER-01** — Cucumber, Non-bitter fruit while retaining bitter-leaf defense
  - symbols that do not resolve: `Bt`, `CsBt`, `Bl`, `CsBl`, `Bi`, `CsBi`
  - confirmed paper: PMID 25430763 — Plant science. Biosynthesis, regulation, and domestication of bitterness in cucumber.
  - needed: the gene symbol used in that paper for Cucumis sativus, or an NCBI/Ensembl gene ID
- **WAT-BRIX-01** — Watermelon, Sugar accumulation (Brix) & flesh color
  - symbols that do not resolve: `ClTST2`, `ClSWEET3`, `ClLCYB`
  - confirmed paper: PMID 29118248 — A Tonoplast Sugar Transporter Underlies a Sugar Accumulation QTL in Watermelon.
  - needed: the gene symbol used in that paper for Citrullus lanatus, or an NCBI/Ensembl gene ID
- **EGG-BROWN-01** — Eggplant, Non-browning cut fruit flesh
  - symbols that do not resolve: `SmPPO4`, `SmPPO5`, `SmPPO6 (PPO`
  - confirmed paper: PMID 33343607 — Simultaneous CRISPR/Cas9 Editing of Three PPO Genes Reduces Fruit Flesh Browning in Solanum mel
  - needed: the gene symbol used in that paper for Solanum melongena, or an NCBI/Ensembl gene ID
- **EGG-PRICKLE-01** — Eggplant, Prickle-free calyx & stem
  - symbols that do not resolve: `SmLOG1`
  - confirmed paper: PMID 38974191 — PE (Prickly Eggplant) encoding a cytokinin-activating enzyme responsible for the formation of p
  - needed: the gene symbol used in that paper for Solanum melongena, or an NCBI/Ensembl gene ID
- **LET-VITC-01** — Lettuce, Elevated Vitamin C (Ascorbate) via uORF editing
  - symbols that do not resolve: `LsGGP2`, `uORF`
  - confirmed paper: PMID 30080209 — Genome editing of upstream open reading frames enables translational control in plants.
  - needed: the gene symbol used in that paper for Lactuca sativa, or an NCBI/Ensembl gene ID
- **BRA-CLUB-01** — Brassica (Broccoli/Cabbage/Cauliflower), Clubroot resistance (Plasmodiophora brassicae)
  - symbols that do not resolve: `CRa`, `CRb`, `Crr1`, `Crr3`
  - confirmed paper: PMID 23054353 — Molecular characterization of the CRa gene conferring clubroot resistance in Brassica rapa.
  - needed: the gene symbol used in that paper for Brassica oleracea / B. rapa, or an NCBI/Ensembl gene ID
- **ONI-TEAR-01** — Onion, Tearless / low-lachrymatory onion & reduced pungency
  - symbols that do not resolve: `AcLFS`, `AcACS1`, `AcACS2`
  - confirmed paper: PMID 18583530 — Silencing onion lachrymatory factor synthase causes a significant change in the sulfur secondar
  - needed: the gene symbol used in that paper for Allium cepa, or an NCBI/Ensembl gene ID

## Priority 2 — genes resolve, but no reference is confirmed (8 entries)

Identity is settled; functional evidence in the crop is not. A DOI/PMID of a knockout, overexpression, allele-swap or field study would move these straight to Level 2.

- **TOM-MILDEW-01** — Tomato, Powdery mildew resistance (Oidium neolycopersici) (`SlMlo1`)
  - draft cited: "Nekrasov et al., 2017, Sci. Rep." — not findable by author+year+gene
  - unconfirmed topic-search hits, treat with suspicion: `SlMlo1` → PMID 35173386
- **TOM-FIRM-01** — Tomato, Fruit firmness & extended shelf life without ripening block (`SlPL, SlPG2a, SlCEL2`)
  - draft cited: "Uluisik et al., 2016, Nat. Biotechnol.; Wang et al., 2019" — not findable by author+year+gene
- **PEP-VIRUS-01** — Pepper, Potyvirus resistance (PVY, TEV, PepMoV — pvr1 / pvr2) (`eIF4E1, eIFiso4E`)
  - draft cited: "Ruffel et al., 2002; Yoon et al., 2020" — not findable by author+year+gene
  - unconfirmed topic-search hits, treat with suspicion: `eIF4E1` → PMID 42374852; `eIFiso4E` → PMID 42239849
- **CUC-SEX-01** — Cucumber, Gynoecy & parthenocarpy (Female flower dominance & seedless fruit) (`CsACS2, CsACS11, CsWIP1, CsARF8`)
  - draft cited: "Hu et al., 2017, Plant J." — not findable by author+year+gene
- **CAR-COLOR-01** — Carrot, Taproot carotenoid & beta-carotene accumulation (`DcRPGE1, DcAPRR2, DcPSY1, DcPSY2, DcLCYB1`)
  - draft cited: "2024 Functional Carrot Root Color Study; Ellison et al., 2018" — not findable by author+year+gene
  - unconfirmed topic-search hits, treat with suspicion: `DcRPGE1` → PMID 38593056; `DcAPRR2` → PMID 38593056; `DcPSY1` → PMID 38593056; `DcPSY2` → PMID 38593056
- **GAP-COMB-01** — Tomato, Combined Heat + Salinity + Drought Tolerance (Mediterranean field conditions) (`SOS1, SOS2, NHX1, DREB2, HsfA2`)
  - draft cited: "Zandalinas et al., 2021; Rivero et al., 2022" — not findable by author+year+gene
- **GAP-REPRO-01** — Tomato, Pollen & ovule thermotolerance + stigma receptivity under heat waves (`HsfA2, HSP70, APX1`)
  - draft cited: "Rieu et al., 2017; Fragkostefanakis et al., 2016" — not findable by author+year+gene
- **GAP-RECOV-01** — Tomato, Post-stress growth & yield recovery kinetics (`CYP707A1`)
  - draft cited: "2025–2026 Plant Resilience Frameworks" — not findable by author+year+gene
  - unconfirmed topic-search hits, treat with suspicion: `CYP707A1` → PMID 41253153

## Priority 3 — neither gene nor reference resolved (2 entries)

- **CUC-MILDEW-01** — Cucumber, Powdery mildew resistance (`CsMLO1, CsMLO8, CsMLO11`)
  - draft cited: "Berg et al., 2015, Front. Plant Sci."
  - needed: both a gene identifier and a functional paper, or the entry should be retired
- **WAT-SEED-01** — Watermelon, Diploid seedless / stenospermocarpic fruit development (`ClAGL11 (plus ClWIP1 for gynoecy)`)
  - draft cited: "Wang et al., 2021; Cucurbit AGL11 studies"
  - needed: both a gene identifier and a functional paper, or the entry should be retired

## Priority 4 — some genes resolve, others do not (5 entries)

- **TOM-DWARF-01** — Tomato, Compact / dwarf architecture for vertical farming & early harvest (3/5 resolved)
  - unresolved: `SlER`, `SlYTH2`
- **TOM-SIZE-01** — Tomato, Fruit size & pericarp thickness (1/3 resolved)
  - unresolved: `SlKIX8`, `SlKIX9`
- **MEL-SHELF-01** — Melon, Extended shelf life & climacteric ripening control (2/3 resolved)
  - unresolved: `CmNAC-NOR`
- **LET-BOLT-01** — Lettuce, Delayed bolting (premature flowering under summer heat) (1/4 resolved)
  - unresolved: `LsCO`, `LsRGL1`, `LsSOC1`
- **BRA-GLUC-01** — Brassica (Broccoli/Cabbage/Cauliflower), Glucoraphanin / Sulforaphane biofortification (1/4 resolved)
  - unresolved: `HAG1`, `MYB29`, `MYB34`

## Commercial status — a question only a person can answer

The draft claimed Level 1 (commercially validated) for 8 entries. No database can confirm that, so the pipeline never derives Level 1: it needs a named cultivar or released product. If you know the product, name it and the entry becomes Level 1:

- **TOM-GABA-01** — Tomato, High GABA accumulation — draft note: "Commercially validated (Sicilian Rouge High GABA). 2025 studies across 3 elite backgrounds confirm genotype-de"
- **TOM-TYLCV-01** — Tomato, TYLCV resistance (ty-5 locus) — draft note: "Commercially deployed via natural ty-5 allele and validated via CRISPR knockout."
- **PEP-PUNG-01** — Pepper, Pungency control (Capsaicinoid presence vs. non-pungent) — draft note: "Commercially established locus; distinct biochemical step from pAMT."
- **PEP-VIRUS-01** — Pepper, Potyvirus resistance (PVY, TEV, PepMoV — pvr1 / pvr2) — draft note: "Commercially validated; base editing avoids pleiotropic defects of full eIF4E + eIFiso4E loss."
- **CUC-BITTER-01** — Cucumber, Non-bitter fruit while retaining bitter-leaf defense — draft note: "Model case of tissue-specific regulation: knocking out Bt removes fruit bitterness without stripping leaf herb"
- **BRA-CLUB-01** — Brassica (Broccoli/Cabbage/Cauliflower), Clubroot resistance (Plasmodiophora brassicae) — draft note: "Commercially critical, but CR loci are pathotype-specific—NOT universal resistance. Single-gene cultivars brea"
- **BRA-GLUC-01** — Brassica (Broccoli/Cabbage/Cauliflower), Glucoraphanin / Sulforaphane biofortification — draft note: "Commercially validated in high-glucoraphanin broccoli hybrids."
- **ONI-TEAR-01** — Onion, Tearless / low-lachrymatory onion & reduced pungency — draft note: "Commercially validated (tearless onions), but NOT a simple single-gene knockout due to multi-copy LFS genes an"

