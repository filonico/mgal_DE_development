#!/bin/bash

source ~/miniforge3/bin/activate phylo_env 

# build hmm profile of Pfam DEAD/HEAD box domain
hmmbuild 00_input/PF00270.alignment.full.hmm 00_input/PF00270.alignment.full.stk

# generate fasta file to align
cat 00_input/vasa_ref.faa 00_input/OG0007* > 06_vasa_ML/vasa_bivalves_withRefOut.faa

# align with clustalo
clustalo -i 06_vasa_ML/vasa_bivalves_withRefOut.faa --hmm-in=00_input/PF00270.alignment.full.hmm --outfmt=fa --threads=10 -o 06_vasa_ML/vasa_bivalves_withRefOut_aligned.faa -v --output-order=tree-order --force

# trim alignment
trimal -in 06_vasa_ML/vasa_bivalves_withRefOut_aligned.faa -out 06_vasa_ML/vasa_bivalves_withRefOut_aligned_trim04.faa -gt 0.4

# manual removal of low-quality sequences
# re-align again
clustalo -i 06_vasa_ML/vasa_bivalves_withRefOut_reduced.faa --hmm-in=00_input/PF00270.alignment.full.hmm --outfmt=fa --threads=10 -o 06_vasa_ML/vasa_bivalves_withRefOut_reduced_aligned.faa -v --output-order=tree-order --force

# re-trim again
trimal -in 06_vasa_ML/vasa_bivalves_withRefOut_reduced_aligned.faa -out 06_vasa_ML/vasa_bivalves_withRefOut_reduced_aligned_trim04.faa -gt 0.4

# build the ML tree
iqtree2 -s 06_vasa_ML/vasa_bivalves_withRefOut_reduced_aligned_trim04.faa -m MFP -bb 1000 -nstop 100 -T AUTO
