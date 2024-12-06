#! /bin/bash
#
#SBATCH --account=lp_edu_comparativegenomics
#SBATCH --job-name=ibp_cagecleaner
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=2 
#SBATCH --mem=4g
#SBATCH --clusters=genius
#SBATCH --time=03:00:00
#SBATCH --get-user-env

# Load conda module:
module load cluster/genius
module load Miniconda3/4.12.0

# Load the conda environment: (assuming the conda environment exists in the user's configuration)
source /user/leuven/350/vsc35068/.bashrc
conda activate /data/leuven/350/vsc35068/miniconda3/envs/cagecleaner

# Run the program:
python3 main.py ../../data/example_runs/cagecat/STS-cagecat/results/N398V589S066P61_binary.csv ../../data/example_runs/cagecat/STS-cagecat/results/N398V589S066P61_summary.txt 99
