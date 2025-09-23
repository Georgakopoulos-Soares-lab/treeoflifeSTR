#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=200:00:00
#SBATCH --partition=sla-prio
#SBATCH --account=izg5139_sc
#SBATCH --mem=30GB
#SBATCH -J Channi
#SBATCH --output=mout_agat_viral/%x_%j.out
#SBATCH --error=merr_agat_viral/%x_%j.err

echo "Working on bucket ${SLURM_ARRAY_TASK_ID}."
micromamba activate agat
python agatify_gff.py $1 --bucket_id ${SLURM_ARRAY_TASK_ID}
echo "Bucket ${SLURM_ARRAY_TASK_ID} has been processed succesfully!"
