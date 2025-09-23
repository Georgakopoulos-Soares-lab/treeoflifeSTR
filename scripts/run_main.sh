#!/bin/bash

SCHEDULE=$1
MODE=$2
TOTAL_BUCKETS=${3:-5}

for bucket in $(seq 0 $((TOTAL_BUCKETS - 1))); do
    echo "Starting training for bucket $bucket"
    bash submit_main.sh $SCHEDULE $MODE $bucket
done
