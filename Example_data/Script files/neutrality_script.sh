#!/bin/bash
# Slurm Script Input Variables
#SBATCH --job-name=n
#SBATCH --chdir=/lustre07/scratch/1001/split_chr/fst_window
#SBATCH --error=n_%A-%a.error
#SBATCH --output=n_%A-%a.out
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=15G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=0-05:30:00
#SBATCH --array=1-5

# Print
echo "=============="
echo "Evolutionary tests" $SLURM_ARRAY_TASK_ID
echo "=============="
echo ""

# Command
start=$SECONDS
module load cuda/11.4
cd /lustre07/scratch/deshan/1001/split_chr/$SLURM_ARRAY_TASK_ID
./CATE -t parameters.json
./CATE -f parameters.json
./CATE -w parameters.json
./CATE -n parameters.json
./CATE -x parameters.json
./CATE -m parameters.json
./CATE -e parameters.json
end=$SECONDS
echo "duration: $((end-start)) seconds."
