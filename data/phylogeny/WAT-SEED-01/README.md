# WAT-SEED-01 — phylogenetic check of Cla97C08G144840 (2026-11-17)

Sequences (each header checked before use):
- Citrullus lanatus Cla97_v1 (Ensembl REST, rest.ensembl.org): Cla97C08G144840.1 (227 aa), Cla97C03G058570.1 (246), Cla97C01G024020.1 (229), Cla97C02G044160.1 (267) — the last three are the within-species paralogues Ensembl Compara reports at Magnoliopsida level.
- Arabidopsis (Ensembl REST): STK AT4G09960.1 (230), SHP1 AT3G58780.4 (281), SHP2 AT2G42830.2 (248), AG AT4G18960.1 (252), SEP3 AT1G24260.2 (251, outgroup).
- NCBI RefSeq: NP_001234194.1 tomato AGL11 (223), NP_001300859.1 TAGL1 (269), NP_001295225.1 tomato AGAMOUS (249), NP_001268118.1 grape AGL11 (223).
- Dropped: cucumber Csa_5G650625 (Compara one-to-one ortholog) returned only a 114 aa fragment (KGN52702).

Method: MAFFT L-INS-i on full length; trimmed to the columns covering Cla97C08G144840 residues 1–175 (MADS, I and K regions); FastTree LG+CAT; 200 column bootstraps; rooted on SEP3.

Files: full.fa, aln.fa (full alignment), aln_mik.fa (trimmed), tree_rooted.nwk, boot.nwk (200 bootstrap trees), compara_orth.json (Ensembl Compara orthologues of Cla97C08G144840).
