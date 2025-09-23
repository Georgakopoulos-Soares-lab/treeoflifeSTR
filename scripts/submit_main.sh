#!/bin/bash

#SBATCH -J NonBDNA
#SBATCH --partition=standard
#SBATCH --mem=45GB
#SBATCH --time=104:00:00
#SBATCH --account=izg5139_cr_default

# JOBDIR="debug_${SLURM_JOB_NAME}"
# mkdir -p "$JOBDIR"

#SBATCH --output=debug_nonbdna_tempura/slurm-%A_%a.out
#SBATCH --error=debug_nonbdna_tempura/slurm-%A_%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user=nmc6088@psu.edu

mkdir -p debug_nonbdna_STR
SCHEDULE=$1
PATTERN=${2:-"STR"}
BID=$3
# BID=${SLURM_ARRAY_TASK_ID}

mkdir -p "extractions_${PATTERN}"
echo "Processing Bucket ${BID} for PATTERN ${PATTERN}."
python main.py --schedule $SCHEDULE \
	       --bucket_id $BID \
	       --pattern $PATTERN \
	       --outdir "extractions_${PATTERN}"
echo "Bucket ${BID} is complete."
