#!/bin/bash

ASSEMBLIES_DIR="assemblies_sample"
ASSEMBLIES_PATH="assemblies_sample.txt"
TOTAL_BUCKETS=${1:-4}

mkdir -p $ASSEMBLIES_DIR
datasets download genome accession --inputfile accessions_sample.txt --filename assemblies_sample.zip
unzip -n assemblies_sample.zip
cp assemblies_sample/**/*.fna $ASSEMBLIES_DIR
find assemblies_sample -type f -name "*fna" > $ASSEMBLIES_PATH
python scheduling.py --accessions_path $ASSEMBLIES_PATH \
                      --total_buckets $TOTAL_BUCKETS

      
