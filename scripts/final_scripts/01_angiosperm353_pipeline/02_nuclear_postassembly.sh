#!/bin/bash
#SBATCH -D ./path2workDir
#SBATCH --job-name=nuclear_postassembly_flax
#SBATCH --get-user-env
#SBATCH -p light
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=logs/nuclear_postassembly_%j.out
#SBATCH --error=logs/nuclear_postassembly_%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=beatrice.landoni@unimi.it
#SBATCH --account=Meristem-Size

# ── User settings ────────────────────────────────────────────────────────────
TARGET="nuclear_targets.fasta"                 # protein target file
NAMELIST="namelist_nuclear.txt" #same as namelist.txt, but with "_nuclear" appended after each sample name
STATS_DIR="stats_nuclear"
# ─────────────────────────────────────────────────────────────────────────────

mkdir -p "$STATS_DIR"

# Activate conda environment
module load anaconda/2025
source /exa/software/Anaconda/2025/.bashanaconda

source activate base
conda activate /home/users/lucio.conti/.conda/envs/hybseq

conda env list | cat > check_conda_env.txt

# Run hybpiper sample checks
echo "[$(date)] Running hybpiper stats..."
hybpiper stats \
    -t_dna "$TARGET" \
    gene "$NAMELIST" \
    --stats_filename "${STATS_DIR}/hybpiper_stats"

mv seq_lengths.tsv "${STATS_DIR}/"
mv gene_read_counts_all* "${STATS_DIR}/"

echo "[$(date)] Running recovery heatmap..."
hybpiper recovery_heatmap "${STATS_DIR}/seq_lengths.tsv"

mv recovery_heatmap* "${STATS_DIR}/"

echo "[$(date)] Running paralog retriever..."
hybpiper paralog_retriever "$NAMELIST" \
    -t_dna "$TARGET" \
    --fasta_dir_all paralogs_all_nuclear \
    --fasta_dir_no_chimeras paralogs_no_chimeras_nuclear

mv paralog_* "${STATS_DIR}/"
mv paralogs_above* "${STATS_DIR}/"

echo "[$(date)] Done. Review heatmap and paralog outputs before submitting 03_nuclear_retrieve.sh."

conda deactivate
