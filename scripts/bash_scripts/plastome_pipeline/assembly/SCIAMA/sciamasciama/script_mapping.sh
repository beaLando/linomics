#
#PBS -l mem=700mb,walltime=72:00:00,nodes=1:ppn=3
#PBS -o log1.txt 
#PBS -e err1.txt 
#PBS -M up869307@myport.ac.uk 
#PBS -m abe
#


fastqc /mnt/lustre/bland/assembly/Linum_RAW/*.fq.gz -o /mnt/lustre/bland/assembly/Linum_RAW -noextract


for f1 in /mnt/c/Users/LandoniB/Plastome/SCIAMA/*F.fq.gz
do
echo "working with file $f1"

dir="/mnt/c/Users/LandoniB/Plastome/SCIAMA/"
f2=${f1%%F.fq.gz}"R.fq.gz"
f1p=${dir}/$(basename -s .fq.gz $f1)_P.filtered.fq.gz
f1u=${dir}/$(basename -s .fq.gz $f1)_U.filtered.fq.gz 
f2p=${dir}/$(basename -s .fq.gz $f2)_P.filtered.fq.gz 
f2u=${dir}/$(basename -s .fq.gz $f2)_U.filtered.fq.gz 

java -jar /home/bland/Programs/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 -trimlog ${dir}/TrimLog $f1 $f2 $f1p $f1u $f2p $f2u ILLUMINACLIP:/mnt/c/Users/LandoniB/Plastome/SCIAMA/refs/adapterIL.fa:2:30:10 AVGQUAL:20 MINLEN:40

done


fastqc /mnt/lustre/bland/assembly/Linum_RAW/*P.filtered.fq.gz -o /mnt/lustre/bland/assembly/Linum_RAW -noextract

multiqc -ip /mnt/lustre/bland/assembly/Linum_RAW/*_fastqc.zip


genome=/mnt/lustre/bland/assembly/Linum_RAW/refs/ref_Lus_NC_036356.fasta
bowtie2-build -f $genome Lus_Chloroplast

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
    bowtie2 -x Lus_Chloroplast -1 $fq1 -2 $fq2 -S $sam
    samtools view -bS $sam > $bam
    samtools sort $bam -o $sorted_bam
    samtools index $sorted_bam
    bcftools mpileup -Ou -f $genome $sorted_bam | bcftools call --ploidy 1 -mv -Oz -o $variants 
    tabix $variants
    rm $fq1 $fq2

    done