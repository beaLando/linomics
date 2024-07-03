#
#PBS -l mem=700mb,walltime=72:00:00,nodes=1:ppn=3
#PBS -o log1.txt 
#PBS -e err1.txt 
#PBS -M up869307@myport.ac.uk 
#PBS -m abe
#


fastqc /mnt/lustre/bland/assembly/Linum_RAW/*.fq.gz -o /mnt/lustre/bland/assembly/Linum_RAW -noextract


for f1 in /mnt/lustre/bland/assembly/Linum_RAW/*F.fq.gz
do
echo "working with file $f1"

dir="/mnt/lustre/bland/assembly/Linum_RAW"
f2=${f1%%F.fq.gz}"R.fq.gz"
f1p=${dir}/$(basename -s .fq.gz $f1)_P.filtered.fq.gz
f1u=${dir}/$(basename -s .fq.gz $f1)_U.filtered.fq.gz 
f2p=${dir}/$(basename -s .fq.gz $f2)_P.filtered.fq.gz 
f2u=${dir}/$(basename -s .fq.gz $f2)_U.filtered.fq.gz 

trimmomatic PE -phred33 -trimlog ${dir}/TrimLog $f1 $f2 $f1p $f1u $f2p $f2u ILLUMINACLIP:${dir}/refs/adapterIL.fa:2:30:10 AVGQUAL:20 MINLEN:40

fastqc /mnt/lustre/bland/assembly/Linum_RAW/*P.filtered.fq.gz -o /mnt/lustre/bland/assembly/Linum_RAW -noextract

done