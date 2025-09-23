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
if [[ -z "$bucket" ]]; then
  echo "Error: SLURM_ARRAY_TASK_ID is not set. Please run this script as part of a SLURM array job."
  bucket=$4
fi

echo "Processing bucket ${bucket}..."
python shuffle.py --schedule $schedule --outdir $outdir --bucket $bucket --level $level
echo "Bucket ${bucket} has been processed succesfully."
