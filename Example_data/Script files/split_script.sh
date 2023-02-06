#!/bin/bash
# Slurm Script Input Variables
#SBATCH --job-name=w
#SBATCH --chdir=/lustre07/scratch/1001/split_chr/fst_window
#SBATCH --error=w_%A-%a.error
#SBATCH --output=w_%A-%a.out
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=15G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=0-05:30:00
#SBATCH --array=1-5

# Print
echo "=============="
echo "Split Chromosome" $SLURM_ARRAY_TASK_ID
echo "=============="
echo ""

# Command
start=$SECONDS
module load cuda/11.4
cd /lustre07/scratch/prometheus_test/full/$SLURM_ARRAY_TASK_ID
./CATE -svcf parameters.json
end=$SECONDS
echo "duration: $((end-start)) seconds."
