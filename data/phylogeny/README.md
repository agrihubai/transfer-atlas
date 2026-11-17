# Phylogenetic checks (2026-11-17)

One folder per gene family. Each holds the input proteins (`*.fa`, all from Ensembl REST, rest.ensembl.org; tomato from the `solanum_lycopersicum_gca000188115v5cm` build), the MAFFT L-INS-i alignment (`*.aln`), the alignment with columns >50% gaps removed (`*.trim.fa`), the FastTree LG+CAT tree rooted on the `_out` outgroup (`*.tree.nwk`), and the splits of 200 column-bootstrap trees (`*.boot.json`).

The Ensembl Compara orthologue calls used alongside the trees are in `compara_*.txt`, `compara_kix.json` and `ER/compara_tomato.txt`. Compara is an independent database answer; the trees are ours. Where both agree, the atlas records the gene; where only Compara speaks (KIX8, KIX9), it says so.

`WAT-SEED-01/` is the earlier AGL11/STK check, with its own README.
