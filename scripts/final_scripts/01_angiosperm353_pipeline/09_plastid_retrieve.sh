#!/bin/bash
#SBATCH -D ./path2workDir
#SBATCH --job-name=plastid_retrieve_flax
#SBATCH --get-user-env
#SBATCH -p light
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=logs/retrieve_%j.out
#SBATCH --error=logs/retrieve_%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=beatrice.landoni@unimi.it
#SBATCH --account=Meristem-Size

# ── User settings ────────────────────────────────────────────────────────────
TARGET="plastid_targets.faa"
NAMELIST="namelist_plastid.txt"
STATS_DIR="./stats_plastid"
OUT_DIR="gene_sequences_plastid25"
# Adjust --filter_by parameters after reviewing recovery_heatmap
FILTER_STAT="GenesAt50pct"
FILTER_OP="greater"
FILTER_VAL=19
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
    -t_aa "$TARGET" \
    --sample_names "$NAMELIST" \
    --stats_file "${STATS_DIR}/hybpiper_stats.tsv" \
    --filter_by "$FILTER_STAT" "$FILTER_OP" "$FILTER_VAL" \
    --fasta_dir "$OUT_DIR"

echo "[$(date)] Retrieval complete. Sequences written to: $OUT_DIR"

conda deactivate
