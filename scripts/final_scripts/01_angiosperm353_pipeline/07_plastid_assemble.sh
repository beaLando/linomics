#!/bin/bash
#SBATCH -D ./path2workDir
#SBATCH --job-name=plastid_assemble_flax
#SBATCH --get-user-env
#SBATCH -p light
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --array=1-28   # one task per sample
#SBATCH --output=logs/assemble_%A_%a.out
#SBATCH --error=logs/assemble_%A_%a.err
#SBATCH --mail-type=end
#SBATCH --mail-user=beatrice.landoni@unimi.it
#SBATCH --account=Meristem-Size

# ── User settings ────────────────────────────────────────────────────────────
TARGET="plastid_targets.faa"                 # protein target file
READS_DIR="./trimmed_fqs"
NAMELIST="namelist.txt"
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p logs

# Activate conda environment
module load anaconda/2025
source /exa/software/Anaconda/2025/.bashanaconda

source activate base
conda activate /home/users/lucio.conti/.conda/envs/hybseq

conda env list | cat > check_conda_env.txt

# ── Per-sample assembly (runs in parallel across array tasks) ─────────────────
NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$NAMELIST")

echo "[$(date)] Starting assembly for sample: $NAME"

hybpiper assemble \
    -t_aa "$TARGET" \
    -r "${READS_DIR}/${NAME}_1_trimmed.fastq.gz" \
        "${READS_DIR}/${NAME}_2_trimmed.fastq.gz" \
    --prefix "$NAME"_plastid \
    --cpu "$SLURM_CPUS_PER_TASK"

echo "[$(date)] Finished assembly for sample: $NAME"

conda deactivate
