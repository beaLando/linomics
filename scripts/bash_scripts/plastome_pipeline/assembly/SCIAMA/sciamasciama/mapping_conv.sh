#
#PBS -l mem=700mb,walltime=72:00:00,nodes=1:ppn=3
#PBS -o logmapping.txt 
#PBS -e errmapping.txt 
#PBS -M up869307@myport.ac.uk 
#PBS -m abe
#

genome=/mnt/lustre/bland/assembly/Linum_RAW/refs/ref_Lus_NC_036356.fasta
bowtie2-build -f $genome /mnt/lustre/bland/assembly/Linum_RAW/refs/Lus_Chloroplast

for fq1 in `ls /mnt/lustre/bland/assembly/Linum_RAW/*F_P.filtered.fq.gz`
    do
    echo "working with file $fq1"

    dir="/mnt/lustre/bland/assembly/Linum_RAW"
    base=$(basename $fq1 "_F_P.filtered.fq.gz")
    echo "base name is $base"

    fq1=${dir}/${base}_F_P.filtered.fq
    fq2=${dir}/${base}_R_P.filtered.fq
    sam=${dir}/${base}.sam
    bam=${dir}/${base}.bam
    sorted_bam=${dir}/${base}.sorted.bam
    variants=${dir}/${base}.calls.vcf.gz

    gunzip -k ${dir}/${base}_F_P.filtered.fq.gz ${dir}/${base}_R_P.filtered.fq.gz
    bowtie2 -x /mnt/lustre/bland/assembly/Linum_RAW/refs/Lus_Chloroplast -1 $fq1 -2 $fq2 -S $sam
    samtools view -bS $sam > $bam
    samtools sort $bam -o $sorted_bam
    samtools index $sorted_bam
    bcftools mpileup -Ou -f $genome $sorted_bam | bcftools call --ploidy 1 -mv -Oz -o $variants 
    tabix $variants
    rm $fq1 $fq2

    done