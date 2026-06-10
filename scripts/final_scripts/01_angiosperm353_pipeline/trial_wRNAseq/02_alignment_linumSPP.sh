# If GDrive does not appear under /mnt/g, you need to mount it with:
# sudo mount -t drvfs G: /mnt/g


# Create alignment subfolder in A353 data directory
mkdir ./data/Angiosperm353/aln_nuclear


# Extract genes for outgroups (euphorbia, phyllanthus) to add to Linum spp.:
## 1. Remove everything except for geneID in fasta headers for both species
## 2. Add species names to each gene header
## 3. Linearize

#cat ./data/Angiosperm353/PAFTOL_raw/oneKP.RHAU.Euphorbia_mesembryanthemifolia.a353.fasta | seqkit replace -p "\s.+"| seqkit grep -nv -f ./data/Angiosperm353/hybpiper/genes2remove.txt | seqkit replace -p $ -r "_euphorbia" | seqkit seq -w 0 | cat > ./data/Angiosperm353/PAFTOL_raw/oneKP.RHAU.Euphorbia_mesembryanthemifolia.a353.clean.fasta
#cat ./data/Angiosperm353/PAFTOL_raw/oneKP.YGAT.Phyllanthus_sp.a353.fasta | seqkit replace -p "\s.+"| seqkit grep -nv -f ./data/Angiosperm353/hybpiper/genes2remove.txt | seqkit replace -p $ -r "_phyllanthus" | seqkit seq -w 0 | cat > ./data/Angiosperm353/PAFTOL_raw/oneKP.YGAT.Phyllanthus_sp.a353.clean.fasta
#cat ./data/Angiosperm353/PAFTOL_raw/oneKP.BVOF.Linum_perenne.a353.fasta | seqkit replace -p "\s.+"| seqkit grep -nv -f ./data/Angiosperm353/hybpiper/genes2remove.txt | seqkit replace -p $ -r "_perenne" | seqkit seq -w 0 | cat > ./data/Angiosperm353/PAFTOL_raw/oneKP.BVOF.Linum_perenne.a353.clean.fasta

#mkdir ./data/Angiosperm353/PAFTOL_raw/outgroups_singlegenes

#cat ./data/Angiosperm353/PAFTOL_raw/oneKP.RHAU.Euphorbia_mesembryanthemifolia.a353.clean.fasta | seqkit split -i -O ./data/Angiosperm353/PAFTOL_raw/outgroups_singlegenes
#cat ./data/Angiosperm353/PAFTOL_raw/oneKP.YGAT.Phyllanthus_sp.a353.clean.fasta | seqkit split -i -O ./data/Angiosperm353/PAFTOL_raw/outgroups_singlegenes#cat ./data/Angiosperm353/PAFTOL_raw/oneKP.BVOF.Linum_perenne.a353.clean.fasta | seqkit split -i -O ./data/Angiosperm353/PAFTOL_raw/outgroups_singlegenes

#for f in `ls ./data/Angiosperm353/hybpiper/gene_sequences/*.FNA`
#do

#gene_name=$(basename -s .FNA $f)
#file_name="./data/Angiosperm353/PAFTOL_raw/outgroups_singlegenes/stdin.part_"$gene_name"_"
#echo ${file_name}

#    cat $file_name*".fasta" >> $f

#done

# Write files for each gene and species: 

##1. linearise fasta files & remove duplicate taxa (extra euphorbia and perenne samples: ERR2040369_nuclear, ERR2040381_nuclear)
##2. cleanse fasta headers from hybpiper notations (these are important for knowing from where hybpiper extracted seqs -e.g. stiched contig or not- but will be troublesome in alignment and tree making)
##3. set aside genes that showed too many paralogs & missing data (you need to go to excel file to create list based on hybpiper stats output)

for f in `ls ./data/Angiosperm353/hybpiper/gene_sequences_nuclear/*.FNA`

do

gene_name=$(basename -s .FNA $f)
echo ${gene_name}

output_name=${gene_name}.fasta

cat $f | seqkit rmdup -n | seqkit replace -p "\s.+" | seqkit seq -w 0 | seqkit grep -v -r -p "ERR2040369|ERR2040381|C95|D04" | cat > "./data/Angiosperm353/aln_nuclear/"$output_name

done


mkdir -p ./data/Angiosperm353/aln_nuclear/not_used

dos2unix ./data/Angiosperm353/hybpiper/genes2remove_nuclear.txt
for file in $(cat ./data/Angiosperm353/hybpiper/genes2remove_nuclear.txt); do mv "$file" ./data/Angiosperm353/aln_nuclear/not_used/; done


# Run mafft alignment for each gene

for raw in `ls ./data/Angiosperm353/aln_nuclear/*.fasta`

do

gene_name=$(basename -s .fasta $raw)
echo ${gene_name}

mafft --thread 5 --preservecase --nuc --localpair --maxiterate 1000 $raw > ./data/Angiosperm353/aln_nuclear/${gene_name}.mft

done


# Trim & linearize clean fastas

for mft in `ls ./data/Angiosperm353/aln_nuclear/*.mft`

do

gene_name=$(basename -s .mft $mft)
echo ${gene_name}

trimal -in $mft -out ./data/Angiosperm353/aln_nuclear/${gene_name}.mft.tal -automated1 -fasta 
cat ./data/Angiosperm353/aln_nuclear/${gene_name}.mft.tal | seqkit seq -w 0 | cat > "./data/Angiosperm353/aln_nuclear/"${gene_name}.tall

done

rm -r ./data/Angiosperm353/aln_nuclear/*.mft.tal


# Create concatenated nexus with partitions based on trimmed alignments AKA supermatrix (download ConcatFasta.py at https://github.com/santiagosnchez/ConcatFasta --can also save to phylyp)
## Output can be used to build phylogenetic trees
python ConcatFasta.py --files ./data/Angiosperm353/aln_nuclear/*.tall --outfile ./data/Angiosperm353/aln_nuclear/supermx353.mft.tall.nexus --part --nexus
