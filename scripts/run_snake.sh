#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=4GB

j=$1

if [[ -d ".snakemake" ]];
then
	echo "SNAKEMAKE <CHANNI>"
	# rm -rf .snakemake
fi

if [[ ! -n ${SSH_CONNECTION} ]];
then
  snakemake --snakefile extract_nonbdna.smk \
                --latency-wait 5 \
                --keep-going \
                --cores $j \
                --keep-incomplete
else
  snakemake --snakefile extract_nonbdna.smk \
	    --keep-incomplete \
	    --keep-going \
	    --latency-wait 45 \
	    --cluster-config cluster_settings.yaml \
	    --cluster "sbatch -p {cluster.partition} \
	    			-t {cluster.time} \
				    --mem={cluster.mem} \
				    -c {cluster.ncpus} \
				    --nodes={cluster.nodes} \
				    -J {cluster.jobName} \
				    -o jobOut/{cluster.jobName}-%j.out \
				    -e jobOut/{cluster.jobName}-%j.err" -j $j
fi
# --rerun-triggers mtime \
