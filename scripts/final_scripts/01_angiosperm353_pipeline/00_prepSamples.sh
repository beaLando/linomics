#!/bin/bash
#SBATCH -D ./path2workDir
#SBATCH --job-name=prep_samples_flax
#SBATCH --get-user-env
#SBATCH -p light
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/addNew_%A_%a.out
#SBATCH --error=logs/addNew_%A_%a.err
#SBATCH --mail-type=end
#SBATCH --mail-user=beatrice.landoni@unimi.it
#SBATCH --account=Meristem-Size

mkdir -p logs


# CONDA
module load anaconda/2025
source /exa/software/Anaconda/2025/.bashanaconda

source activate base
conda activate /home/users/lucio.conti/.conda/envs/gh_variants

conda env list | cat > check_conda_env.txt


# GET DATA from PAFTOL
## NB: SraRunTable_downloadLinks.txt can be created from excel file SraRunTable.xlsx in same scripts folder (Table 2)

cd ./angiosperm353_paftol_fqs/ #I prefer to be inside directory to avoid problems with wget and parallel command and paths

wget -i SraRunTable_downloadLinks.txt

nohup parallel --joblog parallel_jobs.log --nice 10 'bzcat {} | gzip -c > {.}.gz' ::: *bz2 & 
rm -r *.bz2

cd ..


# QC and TRIMMING PAFTOL & NEWLY SEQUENCED SAMPLES
fastqc ./angiosperm353_paftol_fqs/*fastq.gz -o ./angiosperm353_paftol_fqs -t 10 -noextract

for first in `ls ./angiosperm353_paftol_fqs/*_1.fastq.gz`;
do
    second=`echo ${first} | sed 's/_1/_2/g'`
    first_out=${first%".fastq.gz"}
    second_out=${second%".fastq.gz"}

    fastp --thread 9 --in1 ${first} --in2 ${second} --out1 ${first_out}_trimmed.fastq.gz --out2 ${second_out}_trimmed.fastq.gz --detect_adapter_for_pe --trim_poly_g --trim_poly_x --html fastp_paftol.html --json fastp_paftol.json
done

mv fastp_paftol* ./angiosperm353_paftol_fqs/


for first in `ls ./angiosperm353_juan_fqs/*_1.fastq.gz`;
do
    second=`echo ${first} | sed 's/_1/_2/g'`
    first_out=${first%".fastq.gz"}
    second_out=${second%".fastq.gz"}

    fastp --thread 9 --in1 ${first} --in2 ${second} --out1 ${first_out}_trimmed.fastq.gz --out2 ${second_out}_trimmed.fastq.gz --detect_adapter_for_pe --trim_poly_g --trim_poly_x --html fastp_juan.html --json fastp_juan.json
done

mv fastp_juan* ./angiosperm353_juan_fqs/

fastqc ./angiosperm353_juan_fqs/*_trimmed.fastq.gz -o ./angiosperm353_juan_fqs -t 9 -noextract 
fastqc ./angiosperm353_paftol_fqs/*_trimmed.fastq.gz -o ./angiosperm353_paftol_fqs -t 9 -noextract

multiqc -ip ./angiosperm353_juan_fqs/*trimmed_fastqc.zip ./angiosperm353_paftol_fqs/*trimmed_fastqc.zip -o ./multiqc_angiosperm353_juan_paftol

conda deactivate

