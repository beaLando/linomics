#!/bin/bash
#SBATCH -D ./path2workDir
#SBATCH --job-name=nuclear_assemble_flax
#SBATCH --get-user-env
#SBATCH -p light
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --array=1-28   # one task per sample
#SBATCH --output=logs/nuclear_assemble_%A_%a.out
#SBATCH --error=logs/nuclear_assemble_%A_%a.err
#SBATCH --mail-type=end
#SBATCH --mail-user=beatrice.landoni@unimi.it
#SBATCH --account=Meristem-Size

# ── User settings ────────────────────────────────────────────────────────────
TARGET="nuclear_targets.fasta"                 # protein target file
READS_DIR="./trimmed_fqs"
NAMELIST="namelist.txt"
# ─────────────────────────────────────────────────────────────────────────────

# Activate conda environment
module load anaconda/2025
source /exa/software/Anaconda/2025/.bashanaconda

source activate base
conda activate /home/users/lucio.conti/.conda/envs/hybseq

conda env list | cat > check_conda_env.txt

# PREPARE FOLDERS
mkdir -p logs
mkdir -p trimmed_fqs

mv ./angiosperm353_juan_fqs/*_trimmed.fastq.gz ./trimmed_fqs
mv ./angiosperm353_paftol_fqs/*_trimmed.fastq.gz ./trimmed_fqs

rename 's/-read_/_/g' ./trimmed_fqs/*

# CREATE LIST OF ACCESSION NAMES
for fq in `ls ./trimmed_fqs/*_1_trimmed.fastq.gz`;
do
    fl=$(basename -- "$fq")
    accession=${fl%"_1_trimmed.fastq.gz"}
    echo $accession | cat >> namelist.txt
done

# GENERATE TARGET GENE FILE
### see here https://hackmd.io/@mossmatters/Sk82Jujhc#15-filtering-and-correcting-your-target-file-with-%60hybpiper-fix_targetfile%60
### and see here: https://github.com/chrisjackson-pellicle/NewTargets > get target file and python script from here
### for my purpose I targeted reference sequences from linaceae, phyllantaceae, and euphorbiaceae

python ./ang353_NewTargets/filter_mega353.py ./ang353_NewTargets/mega353.fasta ./ang353_NewTargets/select_file.txt
mv ./ang353_NewTargets/filtered_target_file.fasta ./nuclear_targets.fasta

dos2unix nuclear_targets.fasta
dos2unix namelist.txt

# RUN HYBPIPER
NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$NAMELIST")

echo "[$(date)] Starting assembly for sample: $NAME"

hybpiper assemble \
    -t_dna "$TARGET" \
    -r "${READS_DIR}/${NAME}_1_trimmed.fastq.gz" \
        "${READS_DIR}/${NAME}_2_trimmed.fastq.gz" \
    --prefix "$NAME"_nuclear \
    --cpu "$SLURM_CPUS_PER_TASK"

echo "[$(date)] Finished assembly for sample: $NAME"

conda deactivate
