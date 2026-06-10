#!/bin/bash
#SBATCH -D ./path2workDir
#SBATCH --job-name=plastid_postassembly_flax
#SBATCH --get-user-env
#SBATCH -p light
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=logs/postassembly_%j.out
#SBATCH --error=logs/postassembly_%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=beatrice.landoni@unimi.it
#SBATCH --account=Meristem-Size

# ── User settings ────────────────────────────────────────────────────────────
TARGET="plastid_targets.faa"                 # protein target file
NAMELIST="namelist_plastid.txt"
STATS_DIR="./stats_plastid"
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p logs "$STATS_DIR"

# Activate conda environment
module load anaconda/2025
source /exa/software/Anaconda/2025/.bashanaconda

source activate base
conda activate /home/users/lucio.conti/.conda/envs/hybseq

conda env list | cat > check_conda_env.txt

# Run hybpiper sample checks
echo "[$(date)] Running hybpiper stats..."
hybpiper stats \
    -t_aa "$TARGET" \
    gene "$NAMELIST" \
    --stats_filename "${STATS_DIR}/hybpiper_stats"

mv seq_lengths.tsv "${STATS_DIR}/"
mv gene_read_counts_all* "${STATS_DIR}/"

echo "[$(date)] Running recovery heatmap..."
hybpiper recovery_heatmap "${STATS_DIR}/seq_lengths.tsv"

mv recovery_heatmap* "${STATS_DIR}/"

echo "[$(date)] Running paralog retriever..."
hybpiper paralog_retriever "$NAMELIST" \
    -t_aa "$TARGET" \
    --fasta_dir_all paralogs_all_plastid \
    --fasta_dir_no_chimeras paralogs_no_chimeras_plastid

mv paralog_* "${STATS_DIR}/"
mv paralogs_above* "${STATS_DIR}/"

echo "[$(date)] Done. Review heatmap and paralog outputs before submitting 09_plastid_retrieve.sh."

conda deactivate
