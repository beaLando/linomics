#
#PBS -l mem=1000mb,walltime=48:00:00,nodes=1:ppn=3
#PBS -o logquality.txt 
#PBS -e errquality.txt 
#PBS -M up869307@myport.ac.uk 
#PBS -m abe
#

module purge
module load system/intel64
module load apps/trimmomatic/0.38

for f1 in `ls /mnt/lustre/bland/assembly/Linum_RAW/*F.fq.gz`
   do
   echo "working with file $f1"

   dir="/mnt/lustre/bland/assembly/Linum_RAW"
   base=$(basename $f1 "_F.fq.gz")
   echo "base name is $base"

   f1=${dir}/${base}_F.fq.gz
   f2=${dir}/${base}_R.fq.gz
   echo "f1 is $f1"
   echo "f2 is $f2"

   f1p=${dir}/${base}_F_P.filtered.fq.gz
   f1u=${dir}/${base}_F_U.filtered.fq.gz
   f2p=${dir}/${base}_R_P.filtered.fq.gz
   f2u=${dir}/${base}_R_U.filtered.fq.gz

   trimmomatic PE -phred33 -trimlog TrimLog $f1 $f2 $f1p $f1u $f2p $f2u ILLUMINACLIP:/mnt/lustre/bland/assembly/Linum_RAW/refs/adapterIL.fa:2:30:10 AVGQUAL:20 MINLEN:40

   done