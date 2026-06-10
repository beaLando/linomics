#!/bin/bash
#SBATCH -D /home/users/lucio.conti/Meristem-Size/bea/flax/Lus_Plastome_paper2024/Lus_Plastome_paper2024
#SBATCH --job-name=gtrg_beast_plastid_data
#SBATCH --get-user-env
#SBATCH -p light
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/gtrgModel5_beast_plastidCorr_data_%A_%a.out
#SBATCH --error=logs/gtrgModel5_beast_plastidCorr_data_%A_%a.err
#SBATCH --mail-type=end
#SBATCH --mail-user=beatrice.landoni@unimi.it
#SBATCH --account=Meristem-Size

mkdir -p logs

# Activate conda environment
module load anaconda/2025
source /exa/software/Anaconda/2025/.bashanaconda

source activate base
conda activate /home/users/lucio.conti/.conda/envs/bestia

conda env list | cat > check_conda_env.txt

#beast -threads 16 beast_analysis/run_empty_plastid.xml
beast -threads 16 beast_analysis/run_plastid.xml

conda deactivate
