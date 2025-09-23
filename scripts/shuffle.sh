#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=20GB
#SBATCH -J Shuffler
#SBATCH --out=shuffle_%x_%j.out
#SBATCH --err=shuffle_%x_%j.err

schedule=${1:-schedule.json}
outdir=${2:-STR_shuffler}
level=${3:-2}

mkdir -p ${outdir}
bucket=${SLURM_ARRAY_TASK_ID}
echo "Processing bucket ${bucket}..."
python shuffle.py --schedule $schedule --outdir $outdir --bucket $bucket --level $level
echo "Bucket ${bucket} has been processed succesfully."
