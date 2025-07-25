# MANUAL STEPS
## 1. Download matk (MK099064), ndhf (MK090451), trnl (MK066884) for L. villarianum from ncbi based on Maguilla paper
## 2. Open plastome alignment and downloaded villarianum sequences in Aliview
## 3. In Aliview align each villarianum sequence to plastome alignment using muscle in Aliview
## 4. Save to fasta aligned regions

# CONCATENATE GENE ALIGNMENTS AND SAVE TO NEXUS (download ConcatFasta.py at https://github.com/santiagosnchez/ConcatFasta)
python ConcatFasta.py --files ./data/Lb_plastomes_raw/Consensus/plastid_genes_extractedVillarianum/*_alignment.rtf.txt --outfile ./data/Lb_plastomes_raw/Consensus/plastid_genes_extractedVillarianum/plastid_genes_mx.nexus --part --nexus

# PHYLOGENETIC TREE
iqtree -T AUTO -s ./data/Lb_plastomes_raw/Consensus/plastid_genes_extractedVillarianum/plastid_genes_mx.nexus -pre ./data/Lb_plastomes_raw/Consensus/plastid_genes_extractedVillarianum/plastid_genes_mx -bb 1000 -wbtl -o "L88"

pxrr -t ./data/Lb_plastomes_raw/Consensus/plastid_genes_extractedVillarianum/plastid_genes_mx.treefile -g L88 -o ./data/Lb_plastomes_raw/Consensus/plastid_genes_extractedVillarianum/plastid_genes_mx.phyxed.treefile