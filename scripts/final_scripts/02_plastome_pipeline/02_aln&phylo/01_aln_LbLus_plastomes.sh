conda activate hybpiper

# Run MAFFT and trimAL
mafft --auto --nuc ./mapping/all/all_plastomes.linear.fasta > ./mapping/all/all_plastomes.aligned.fasta ##because sequences are long, I cannot do local alignment as for Angiosperm353

trimal -in ./mapping/all/all_plastomes.aligned.fasta -out ./mapping/all/all_plastomes.aligned.trimmed.fasta -automated1 -fasta 


# IQtree
iqtree -T AUTO -s ./mapping/all/all_plastomes.aligned.trimmed.fasta -m GTR+I+G -bb 1000 -wbtl -o "L88"

pxrr -t ./data/Lb_plastomes_raw/Consensus/aln/all/all_plastomes.aligned.trimmed.fasta.contree -g L88 -o ./data/Lb_plastomes_raw/Consensus/aln/all/all_plastomes.aligned.trimmed.fasta.phyxed.contree



# REROOT also NUCLEAR PHYLOGENY YANN
pxrr -t ./data/Lb_plastomes_raw/Consensus/genetic_distances_tree_nuclear.nwk -g L88 -o ./data/Lb_plastomes_raw/Consensus/genetic_distances_tree_nuclear.phyxed.nwk


# RUN BEAST on CHLOROPLAST PHYLOGENY (point & click)
