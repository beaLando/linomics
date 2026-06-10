# ALIGNMENT, TREEs, AND TREEPL DATING of PLASTID-gene SUPERMATRIX PHYLOGENY
## NB: astral is not considered for plastid genes, since they should be considered as bits from one single large locus


# CREATE ALIGNMENT & TREE SUBFOLDERS
mkdir -p ./data/Angiosperm353/aln_plastid/
#mkdir -p ./data/Angiosperm353/aln_plastid25/ #this recovery % was not used in the end
mkdir -p ./data/Angiosperm353/iqtrees/plastid
#mkdir -p ./data/Angiosperm353/iqtrees/plastid25

conda activate phylo


# WRITE FILES FOR EACH GENE AND SPECIES
## 1. linearise fasta files 
## 2. cleanse fasta headers from hybpiper notations (these are important for knowing from where hybpiper extracted seqs -e.g. stiched contig or not- but will be troublesome in alignment and tree making)
## 3. set aside genes that showed too many paralogs & missing data (you need to go to excel file to create list based on hybpiper stats output)

## 50 PERC. SEQ RECOVERY
for f in `ls ./data/Angiosperm353/hybpiper/gene_sequences_plastid/*.FNA`

do

gene_name=$(basename -s .FNA $f)
echo ${gene_name}

output_name=${gene_name}.fasta

cat $f | seqkit rmdup -n | seqkit replace -p "\s.+" | seqkit seq -w 0 | seqkit grep -v -r -p "D06|D04|D08|D03|D02" | cat > "./data/Angiosperm353/aln_plastid/"$output_name

done


mkdir ./data/Angiosperm353/aln_plastid/not_used

sudo dos2unix ./data/Angiosperm353/hybpiper/genes2remove_plastid.txt
for file in $(cat ./data/Angiosperm353/hybpiper/genes2remove_plastid.txt); do mv "$file" ./data/Angiosperm353/aln_plastid/not_used/; done


## 25 PERC. SEQ RECOVERY
#for f in `ls ./data/Angiosperm353/hybpiper/gene_sequences_plastid25/*.FNA`

#do

#gene_name=$(basename -s .FNA $f)
#echo ${gene_name}

#output_name=${gene_name}.fasta

#cat $f | seqkit rmdup -n | seqkit replace -p "\s.+" | seqkit seq -w 0 | cat > "./data/Angiosperm353/aln_plastid25/"$output_name

#done


#mkdir ./data/Angiosperm353/aln_plastid25/not_used

#sudo dos2unix ./data/Angiosperm353/hybpiper/genes2remove_plastid25.txt
#for file in $(cat ./data/Angiosperm353/hybpiper/genes2remove_plastid25.txt); do mv "$file" ./data/Angiosperm353/aln_plastid25/not_used/; done



# RUN MAFFT for EACH GENE
## 50 PERC. SEQ RECOVERY
for raw in `ls ./data/Angiosperm353/aln_plastid/*.fasta`

do

gene_name=$(basename -s .fasta $raw)
echo ${gene_name}

mafft --thread 5 --preservecase --nuc --localpair --maxiterate 1000 $raw > ./data/Angiosperm353/aln_plastid/${gene_name}.mft

done

## 25 PERC. SEQ RECOVERY
#for raw in `ls ./data/Angiosperm353/aln_plastid25/*.fasta`

#do

#gene_name=$(basename -s .fasta $raw)
#echo ${gene_name}

#mafft --thread 5 --preservecase --nuc --localpair --maxiterate 1000 $raw > ./data/Angiosperm353/aln_plastid25/${gene_name}.mft

#done



# TRIM & LINEARIZE CLEAN ALIGNMENTs
## 50 PERC. SEQ RECOVERY
for mft in `ls ./data/Angiosperm353/aln_plastid/*.mft`

do

gene_name=$(basename -s .mft $mft)
echo ${gene_name}

trimal -in $mft -out ./data/Angiosperm353/aln_plastid/${gene_name}.mft.tal -automated1 -fasta 
cat ./data/Angiosperm353/aln_plastid/${gene_name}.mft.tal | seqkit seq -w 0 | cat > "./data/Angiosperm353/aln_plastid/"${gene_name}.tall

done

rm -r ./data/Angiosperm353/aln_plastid/*.mft.tal


## 25 PERC. SEQ RECOVERY
#for mft in `ls ./data/Angiosperm353/aln_plastid25/*.mft`

#do

#gene_name=$(basename -s .mft $mft)
#echo ${gene_name}

#trimal -in $mft -out ./data/Angiosperm353/aln_plastid25/${gene_name}.mft.tal -automated1 -fasta 
#cat ./data/Angiosperm353/aln_plastid25/${gene_name}.mft.tal | seqkit seq -w 0 | cat > "./data/Angiosperm353/aln_plastid25/"${gene_name}.tall

#done

#rm -r ./data/Angiosperm353/aln_plastid25/*.mft.tal



# CREATE CONCATENATED MATRICES (plastome is one big locus)
## Creates concatenated nexus with partitions based on trimmed alignments AKA supermatrix (download ConcatFasta.py at https://github.com/santiagosnchez/ConcatFasta --can also save to phylyp)
## Output can be used to build phylogenetic trees
python ConcatFasta.py --files ./data/Angiosperm353/aln_plastid/*.tall --outfile ./data/Angiosperm353/aln_plastid/supermx353.mft.tall.nexus --part --nexus
#python ConcatFasta.py --files ./data/Angiosperm353/aln_plastid25/*.tall --outfile ./data/Angiosperm353/aln_plastid25/supermx353.mft.tall.nexus --part --nexus



# MAKE SUPERMATRIX PLASTID TREEs
iqtree -nt 5 -s ./data/Angiosperm353/aln_plastid/supermx353.mft.tall.nexus -m GTR+G -bb 1000 -wbtl -pre ./data/Angiosperm353/iqtrees/plastid/supermx353 
#iqtree -nt 5 -s ./data/Angiosperm353/aln_plastid25/supermx353.mft.tall.nexus -m GTR+G -bb 1000 -wbtl -pre ./data/Angiosperm353/iqtrees/plastid25/supermx353 

### ROOT IT
pxrr -t ./data/Angiosperm353/iqtrees/plastid/supermx353.treefile -g ERR2040368_plastid -o ./data/Angiosperm353/iqtrees/plastid/supermx353.phyxed.treefile #ERR2040368 is euphorbia
#pxrr -t ./data/Angiosperm353/iqtrees/plastid25/supermx353.treefile -g ERR2040368_plastid -o ./data/Angiosperm353/iqtrees/plastid25/supermx353.phyxed.treefile #ERR2040368 is euphorbia



# DATE TREEs
## PRIMING
treePL ./data/Angiosperm353/iqtrees/plastid/config_supermx_s0_linumCN.txt > ./data/Angiosperm353/iqtrees/plastid/config_supermx_s0_linumCN.log

treePL ./data/Angiosperm353/iqtrees/plastid/config_supermx_s100_linumCN.txt > ./data/Angiosperm353/iqtrees/plastid/config_supermx_s100_linumCN.log

## DATING
treePL ./data/Angiosperm353/iqtrees/plastid/config_supermx_s0_linumCN_primed.txt > ./data/Angiosperm353/iqtrees/plastid/config_supermx_s0_linumCN_primed.log

treePL ./data/Angiosperm353/iqtrees/plastid/config_supermx_s100_linumCN_primed.txt > ./data/Angiosperm353/iqtrees/plastid/config_supermx_s100_linumCN_primed.log

