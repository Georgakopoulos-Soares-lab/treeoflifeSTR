#!/bin/bash

ASSEMBLIES_DIR="assemblies_sample"
ASSEMBLIES_PATH="assemblies_sample.txt"
TOTAL_BUCKETS=${1:-4}
mkdir -p $ASSEMBLIES_DIR

datasets download genome accession --inputfile accessions_sample.txt 
unzip -n ncbi_dataset.zip
mv ncbi_dataset/data/**/*.fna $ASSEMBLIES_DIR
find $ASSEMBLIES_DIR -type f -name "*fna" > $ASSEMBLIES_PATH

if [[ -f $ASSEMBLIES_PATH ]]; 
then
	python scheduling.py --accessions_path $ASSEMBLIES_PATH \
                      --total_buckets $TOTAL_BUCKETS
	echo "Success!"
else
	echo "An error has occurred."
fi

      
