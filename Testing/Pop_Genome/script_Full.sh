#!/bin/bash
# Slurm Script Input Variables
#SBATCH --job-name=popGenome
#SBATCH --chdir=/CATE/Testing/Pop_Genome/R_script/
#SBATCH --error=popGenome.error
#SBATCH --output=popGenome.out
#SBATCH --mem=15G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --partition=test

# Print
echo "=============="
echo "Running Running test framework for PopGenome"
echo "=============="
echo ""

# Command
start=$SECONDS
module load R/3.6.2
Rscript script_sample.R
end=$SECONDS
echo "duration: $((end-start)) seconds."
