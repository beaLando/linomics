#!/bin/bash
#SBATCH -D /home/users/lucio.conti/Meristem-Size/bea/flax/Lus_Plastome_paper2024/Lus_Plastome_paper2024
#SBATCH --job-name=gtrgModel_plastomes_flax
#SBATCH --get-user-env
#SBATCH -p light
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/gtrgModel_plastomes_%A_%a.out
#SBATCH --error=logs/gtrgModel_plastomes_%A_%a.err
#SBATCH --mail-type=end
#SBATCH --mail-user=beatrice.landoni@unimi.it
#SBATCH --account=Meristem-Size

mkdir -p logs

# Activate conda environment
module load anaconda/2025
source /exa/software/Anaconda/2025/.bashanaconda

source activate base
conda activate /home/users/lucio.conti/.conda/envs/hybseq

conda env list | cat > check_conda_env.txt

# MAFFT & TRIMal
mafft --auto --nuc ./mapping/all/all_plastomes.linear.fasta > ./mapping/all/all_plastomes.aligned.fasta ##because sequences are long, I cannot do local alignment as for Angiosperm353

trimal -in ./mapping/all/all_plastomes.aligned.fasta -out ./mapping/all/all_plastomes.aligned.trimmed.fasta -automated1 -fasta 


# IQTREE & ROOTING
iqtree -T AUTO -s ./mapping/all/all_plastomes.aligned.trimmed.fasta -m GTR+G -bb 1000 -wbtl -o "L88" -pre ./mapping/all/all_plastomes_gtrg.aligned.trimmed
pxrr -t ./mapping/all/all_plastomes_gtrg.aligned.trimmed.treefile -g L88 -o ./mapping/all/all_plastomes_gtrg.aligned.trimmed.phyxed.treefile


# ALSO ROOT NUCLEAR GENETIC DISTANCE TREE BASED ON YANN'S ANALYSIS
pxrr -t ./data/Lb_plastomes_raw/Consensus/genetic_distances_tree_nuclear.nwk -g L88 -o ./data/Lb_plastomes_raw/Consensus/genetic_distances_tree_nuclear.phyxed.nwk

conda deactivate