#!/bin/bash
#SBATCH -D ./path2workDir
#SBATCH --job-name=nuclear_retrieve_flax
#SBATCH --get-user-env
#SBATCH -p light
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=logs/nuclear_retrieve_%j.out
#SBATCH --error=logs/nuclear_retrieve_%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=beatrice.landoni@unimi.it
#SBATCH --account=Meristem-Size

# ── User settings ────────────────────────────────────────────────────────────
TARGET="nuclear_targets.fasta"
NAMELIST="namelist_nuclear.txt"
STATS_DIR="stats_nuclear"
OUT_DIR="gene_sequences_nuclear"
# Adjust --filter_by parameters after reviewing recovery_heatmap and
# paralog_retriever outputs from script 02
FILTER_STAT="GenesAt50pct"
FILTER_OP="greater"
FILTER_VAL=49
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p logs "$OUT_DIR"

# Activate conda environment
module load anaconda/2025
source /exa/software/Anaconda/2025/.bashanaconda

source activate base
conda activate /home/users/lucio.conti/.conda/envs/hybseq

conda env list | cat > check_conda_env.txt

# Retrieve Sequences
echo "[$(date)] Retrieving sequences..."
hybpiper retrieve_sequences dna \
    -t_dna "$TARGET" \
    --sample_names "$NAMELIST" \
    --stats_file "${STATS_DIR}/hybpiper_stats.tsv" \
    --filter_by "$FILTER_STAT" "$FILTER_OP" "$FILTER_VAL" \
    --fasta_dir "$OUT_DIR"

echo "[$(date)] Retrieval complete. Sequences written to: $OUT_DIR"

conda deactivate
