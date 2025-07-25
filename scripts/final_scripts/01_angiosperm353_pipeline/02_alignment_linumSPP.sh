# If GDrive does not appear under /mnt/g, you need to mount it with:
# sudo mount -t drvfs G: /mnt/g


# Create alignment subfolder in A353 data directory
#mkdir ./data/Angiosperm353/aln/


# Extract genes for outgroups (euphorbia, phyllanthus) to add to Linum spp.:
## 1. Remove everything except for geneID in fasta headers for both species
## 2. Add species names to each gene header
## 3. Linearize

cat ./data/Angiosperm353/PAFTOL_raw/oneKP.RHAU.Euphorbia_mesembryanthemifolia.a353.fasta | seqkit replace -p "\s.+"| seqkit grep -nv -f ./data/Angiosperm353/hybpiper/genes2remove.txt | seqkit replace -p $ -r "_euphorbia" | seqkit seq -w 0 | cat > ./data/Angiosperm353/PAFTOL_raw/oneKP.RHAU.Euphorbia_mesembryanthemifolia.a353.clean.fasta
cat ./data/Angiosperm353/PAFTOL_raw/oneKP.YGAT.Phyllanthus_sp.a353.fasta | seqkit replace -p "\s.+"| seqkit grep -nv -f ./data/Angiosperm353/hybpiper/genes2remove.txt | seqkit replace -p $ -r "_phyllanthus" | seqkit seq -w 0 | cat > ./data/Angiosperm353/PAFTOL_raw/oneKP.YGAT.Phyllanthus_sp.a353.clean.fasta
cat ./data/Angiosperm353/PAFTOL_raw/oneKP.BVOF.Linum_perenne.a353.fasta | seqkit replace -p "\s.+"| seqkit grep -nv -f ./data/Angiosperm353/hybpiper/genes2remove.txt | seqkit replace -p $ -r "_perenne" | seqkit seq -w 0 | cat > ./data/Angiosperm353/PAFTOL_raw/oneKP.BVOF.Linum_perenne.a353.clean.fasta

#mkdir ./data/Angiosperm353/PAFTOL_raw/outgroups_singlegenes

cat ./data/Angiosperm353/PAFTOL_raw/oneKP.RHAU.Euphorbia_mesembryanthemifolia.a353.clean.fasta | seqkit split -i -O ./data/Angiosperm353/PAFTOL_raw/outgroups_singlegenes
cat ./data/Angiosperm353/PAFTOL_raw/oneKP.YGAT.Phyllanthus_sp.a353.clean.fasta | seqkit split -i -O ./data/Angiosperm353/PAFTOL_raw/outgroups_singlegenescat ./data/Angiosperm353/PAFTOL_raw/oneKP.BVOF.Linum_perenne.a353.clean.fasta | seqkit split -i -O ./data/Angiosperm353/PAFTOL_raw/outgroups_singlegenes

for f in `ls ./data/Angiosperm353/hybpiper/gene_sequences/*.FNA`
do

gene_name=$(basename -s .FNA $f)
file_name="./data/Angiosperm353/PAFTOL_raw/outgroups_singlegenes/stdin.part_"$gene_name"_"
echo ${file_name}

    cat $file_name*".fasta" >> $f

done

# Write files for each gene and species: 

##1. linearise fasta files 
##2. cleanse fasta headers from hybpiper notations (these are important for knowing from where hybpiper extracted seqs -e.g. stiched contig or not- but will be troublesome in alignment and tree making)
##3. set aside genes that showed too many paralogs & missing data (you need to go to excel file to create list based on hybpiper stats output)

for f in `ls ./data/Angiosperm353/hybpiper/gene_sequences/*.FNA`

do

gene_name=$(basename -s .FNA $f)
echo ${gene_name}

output_name=${gene_name}.fasta

cat $f | seqkit rmdup -n | seqkit replace -p "\s.+" | seqkit seq -w 0 | cat > "./data/Angiosperm353/aln/"$output_name

done


#mkdir ./data/Angiosperm353/aln/not_used

dos2unix ./data/Angiosperm353/hybpiper/genes2remove.txt
for file in $(cat ./data/Angiosperm353/hybpiper/genes2remove.txt); do mv "$file" ./data/Angiosperm353/aln/not_used/; done


# Run mafft alignment for each gene

for raw in `ls ./data/Angiosperm353/aln/*.fasta`

do

gene_name=$(basename -s .fasta $raw)
echo ${gene_name}

mafft --thread 5 --preservecase --nuc --localpair --maxiterate 1000 $raw > ./data/Angiosperm353/aln/${gene_name}.mft

done


# Trim & linearize clean fastas

for mft in `ls ./data/Angiosperm353/aln/*.mft`

do

gene_name=$(basename -s .mft $mft)
echo ${gene_name}

trimal -in $mft -out ./data/Angiosperm353/aln/${gene_name}.mft.tal -automated1 -fasta 
cat ./data/Angiosperm353/aln/${gene_name}.mft.tal | seqkit replace -p "^.+_" | seqkit seq -w 0 | cat > "./data/Angiosperm353/aln/"${gene_name}.tall

done

rm -r ./data/Angiosperm353/aln/*.mft.tal


# Create concatenated nexus with partitions based on trimmed alignments AKA supermatrix (download ConcatFasta.py at https://github.com/santiagosnchez/ConcatFasta --can also save to phylyp)
## Output can be used to build phylogenetic trees
python ConcatFasta.py --files ./data/Angiosperm353/aln/*.tall --outfile ./data/Angiosperm353/aln/supermx353.mft.tall.nexus --part --nexus

### NB: Previously I had used this method to create concatenated supermatrix, but does not seem good: https://evomics.org/wp-content/uploads/2019/01/Workshop_practical_concatenation_model_testing.pdf (then refine partition file in r)
#java -jar /home/bland/Programs/phyutility/phyutility.jar -concat -in ./data/Angiosperm353/aln/*.mft.tall -out ./data/Angiosperm353/aln/supermx353.mft.tall.nexus
#head -n 3 ./data/Angiosperm353/aln/supermx353.mft.tall.nexus | tail -n 1 | sed 's/^.*\[//g' | sed 's/\s\].*$//g' | sed 's/.mft.tall_gene...\s/=/g'| sed 's/.mft.tall_gene..\s/=/g'| sed 's/.mft.tall_gene.\s/=/g' | tr " " "\n" | sed 's/^/AUTO,\s/g' > ./data/Angiosperm353/aln/supermx353_partitions.txt
#### Turn ambiguity codes into missing in nexus2 (n = ?) ...this works, but I am not sure what is the impact of coding ambiguous site as missing
#sed -e 's/HNCF/HXCF/g:s/phyllanthus/phyllaxthus/g' ./data/Angiosperm353/aln/supermx353.mft.tall.nexus | sed -e '7,${s/n/?/g}' | sed -e '7,${s/x/n/g}' | sed -e 's/HXCF/HNCF/g:s/phyllaxthus/phyllanthus/g' | cat > ./data/Angiosperm353/aln/supermx353.final.mft.tall.nexus
