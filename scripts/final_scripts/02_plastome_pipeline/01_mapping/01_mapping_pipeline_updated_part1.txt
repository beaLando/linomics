# CREATE NECESSARY DIRECTORIES
mkdir ./refs ### LOAD Lus chloroplast genome and illumina adapters to this folder
mkdir ./mapping

# QUALITY CHECK & CLEAN READS
conda activate readNmap

rename 's/Copy of //g' ./raw_fasq_files_plastome_linum/*.gz

fastqc ./raw_fasq_files_plastome_linum/*.fq.gz -o ./raw_fasq_files_plastome_linum -t 6 -noextract

for f1 in ./raw_fasq_files_plastome_linum/*F.fq.gz
do
echo "working with file $f1"

dir="./mapping"
f2=${f1%%F.fq.gz}"R.fq.gz"
logf=${dir}/$(basename -s .fq.gz $f1)
f1p=${dir}/$(basename -s .fq.gz $f1)_P.filtered.fq.gz
f1u=${dir}/$(basename -s .fq.gz $f1)_U.filtered.fq.gz 
f2p=${dir}/$(basename -s .fq.gz $f2)_P.filtered.fq.gz 
f2u=${dir}/$(basename -s .fq.gz $f2)_U.filtered.fq.gz 

fastp -i $f1 -I $f2 -q 20 -u 40 -l 40 --detect_adapter_for_pe --dedup --thread 6 -o $f1p -O $f2p --unpaired1 $f1u --unpaired2 $f2u -j ${logf}.json -h ${logf}.html

done

fastqc ./mapping/*P.filtered.fq.gz -o ./mapping -t 6 -noextract
multiqc -ip -f ./raw_fasq_files_plastome_linum/*_fastqc.zip ./mapping/*_fastqc.zip


# MAP TO REFERENCE
DISPLAY=:0 #for no problem with qualimap
genome=./refs/NC_036356_1_Lus_chloroplast_ref.fasta
bowtie2-build -f $genome Lus_Chloroplast

for fqgz in `ls ./mapping/*F_P.filtered.fq.gz` #around 1min per sample 
    do

    dir="./mapping"
    base=$(basename $fqgz "_F_P.filtered.fq.gz")
    echo "working with sample $base"

    fq1=${dir}/${base}_F_P.filtered.fq
    fq2=${dir}/${base}_R_P.filtered.fq
    sam=${dir}/${base}.sam
    bam=${dir}/${base}.bam
    sorted_bam=${dir}/${base}.sorted.bam
    bam_stats=${dir}/${base}_qualimap
    #variants=${dir}/${base}.calls.vcf.gz

    gunzip -k ${dir}/${base}_F_P.filtered.fq.gz ${dir}/${base}_R_P.filtered.fq.gz
    bowtie2 -x Lus_Chloroplast -1 $fq1 -2 $fq2 -S $sam --threads 3
    samtools view -f0x2 -bS $sam > $bam #f0x2 option is to retain only reads where both mates F and R were kept
    samtools sort $bam -o $sorted_bam
    samtools index $sorted_bam
    qualimap bamqc -bam $sorted_bam -c -nt 6 -outdir $bam_stats

    rm $fq1
    rm $fq2

    done

multiqc -ip -f ./raw_fasq_files_plastome_linum/*_fastqc.zip ./mapping/*_fastqc.zip  ./mapping/*_qualimap


# VARIANT CALLING
for fqgz in `ls ./mapping/*F_P.filtered.fq.gz` 
    do

    dir="./mapping"
    base=$(basename $fqgz "_F_P.filtered.fq.gz")
    echo "working with sample $base"

    sorted_bam=${dir}/${base}.sorted.bam
    variants=${dir}/${base}.calls.vcf.gz
    st=${dir}/${base}.vcf.stats

    bcftools mpileup -Ou -d 450 -f $genome $sorted_bam | bcftools call --ploidy 1 -mv -Oz -o $variants #option d can be modified depending on average coverage depth for example
    tabix $variants

    bcftools stats -F $genome -s - $variants > $st

    done

multiqc -ip -f ./raw_fasq_files_plastome_linum/*_fastqc.zip ./mapping/*_fastqc.zip  ./mapping/*_qualimap ./mapping/*.vcf.stats

# ADD SAMPLE NAME to VCF
for f in `ls ./mapping/*.sorted.bam`
    do
    echo "working with file $f"

    dir="./mapping"
    base=$(basename $f ".sorted.bam")
    out=${dir}/${base}.nhrg.bam
    lb="unknown"
    pl="Illumina"
    pu=${base}_${lb}
    echo -e "@RG\tSM:$base\tPL:$pl\tLB:$lb\tPU:$pu"

    picard AddOrReplaceReadGroups -I $f -O $out -RGLB $lb -RGPL $pl -RGPU $pu -RGSM $base

done

rm -r ./mapping/*.sorted.bam
rm -r ./mapping/*.sorted.bam.bai

for bam in `ls ./mapping/*.nhrg.bam`
    do
    echo "working with file $bam"

    dir="./mapping"
    base=$(basename $bam ".nhrg.bam")
    sorted_bam=${dir}/${base}.nhrg.sorted.bam

    samtools sort $bam -o $sorted_bam
    samtools index $sorted_bam

done

rm -r ./mapping/*.nhrg.bam


# MERGE VCF file (first exclude L84)
mkdir ./mapping/not_used
mkdir ./mapping/all
mv -r ./mapping/L84* ./mapping/not_used/

#bcftools merge -Oz ./mapping/*.calls.vcf.gz > ./mapping/all/all_plastomes.vcf.gz

bcftools mpileup -Ou -d 450 -f $genome ./mapping/*.nhrg.sorted.bam --threads 5 | bcftools call --ploidy 1 -mv -Oz -o ./mapping/all/all_plastomes.vcf.gz #-a AD,DP,SP #-a GQ,GP 
tabix ./mapping/all/all_plastomes.vcf.gz
bcftools stats ./mapping/all/all_plastomes.vcf.gz

vcftools --gzvcf ./mapping/all/all_plastomes.vcf.gz --freq2 --out ./mapping/all/all_plastomes --max-alleles 2
vcftools --gzvcf ./mapping/all/all_plastomes.vcf.gz --depth --out ./mapping/all/all_plastomes
vcftools --gzvcf ./mapping/all/all_plastomes.vcf.gz --site-mean-depth --out ./mapping/all/all_plastomes
vcftools --gzvcf ./mapping/all/all_plastomes.vcf.gz --site-quality --out ./mapping/all/all_plastomes
vcftools --gzvcf ./mapping/all/all_plastomes.vcf.gz --missing-indv --out ./mapping/all/all_plastomes
vcftools --gzvcf ./mapping/all/all_plastomes.vcf.gz --missing-site --out ./mapping/all/all_plastomes
##NB plot with vcf_stats4filtering.R and choose filtering parameters based on that

conda deactivate