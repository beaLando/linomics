# VARIANT FILTERING & CONSENSUS

conda activate readNmap

bcftools filter -e 'MQ<30 || QUAL<30 || DP<10000 || DP>40000' ./mapping/all/all_plastomes.vcf.gz -Oz -o ./mapping/all/all_plastomes.filtered.vcf.gz

tabix ./mapping/all/all_plastomes.filtered.vcf.gz

vcftools --gzvcf ./mapping/all/all_plastomes.filtered.vcf.gz --freq2 --out ./mapping/all/all_plastomes_filtered --max-alleles 2
vcftools --gzvcf ./mapping/all/all_plastomes.filtered.vcf.gz --depth --out ./mapping/all/all_plastomes_filtered
vcftools --gzvcf ./mapping/all/all_plastomes.filtered.vcf.gz --site-mean-depth --out ./mapping/all/all_plastomes_filtered
vcftools --gzvcf ./mapping/all/all_plastomes.filtered.vcf.gz --site-quality --out ./mapping/all/all_plastomes_filtered
vcftools --gzvcf ./mapping/all/all_plastomes.filtered.vcf.gz --missing-indv --out ./mapping/all/all_plastomes_filtered
vcftools --gzvcf ./mapping/all/all_plastomes.filtered.vcf.gz --missing-site --out ./mapping/all/all_plastomes_filtered
##NB plot with vcf_stats4filtering.R and choose filtering parameters based on that

bcftools stats -F $genome -s - ./mapping/all/all_plastomes.filtered.vcf.gz > ./mapping/all/all_plastomes.filtered.vcf.stats

#cat $genome | bcftools consensus ./mapping/all/all_plastomes.filtered.vcf.gz > ./mapping/all/all_plastomes.filtered.consensus.fasta #1959 variants applied


# CONSENSUS

for sample in `bcftools query -l ./mapping/all/all_plastomes.filtered.vcf.gz`
  do

    dir="./mapping"
    base=$(basename $sample ".sorted.bam")
    vcf_filt_ind=${dir}/${base}.vcf.gz
    echo "base name is $base"

    bcftools view -Oz -s $sample -o $vcf_filt_ind ./mapping/all/all_plastomes.filtered.vcf.gz
  
  done

for sample in `bcftools query -l ./mapping/all/all_plastomes.filtered.vcf.gz`
  do

    dir="./mapping"
    base=$(basename $sample ".sorted.bam")
    vcf_filt_ind=${dir}/${base}.vcf.gz
    cons=${dir}/${base}.consensus.fasta
    echo "base name is $base"

    tabix $vcf_filt_ind
    cat $genome | bcftools consensus $vcf_filt_ind > $cons #1958 variants applied
  
  done

conda deactivate


# CONCATENATE FASTAS
conda activate readNmap

grep -e ">" ./mapping/*consensus.fasta #need to rename fasta headers with file names

awk '/^>/{print ">" substr(FILENAME,1,length(FILENAME)-16); next} 1' ./mapping/*consensus.fasta > ./mapping/all/all_plastomes.fasta

cat ./mapping/all/all_plastomes.fasta | seqkit replace -p "\.\/mapping\/" -r "" > ./mapping/all/all_plastomes_rn.fasta

grep -e ">" ./mapping/all/all_plastomes_rn.fasta

cat ./refs/NC_036356_1_Lus_chloroplast_ref.fasta ./mapping/all/all_plastomes_rn.fasta > ./mapping/all/all_plastomes.fasta

cat ./mapping/all/all_plastomes.fasta | seqkit seq -w 0 | cat > ./mapping/all/all_plastomes.linear.fasta

conda deactivate