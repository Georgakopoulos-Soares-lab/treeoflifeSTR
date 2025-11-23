#!/bin/bash

TOTAL_BUCKETS=${1:-2}
SCHEDULE=${2:-schedule.json}
MODE=${3:-STR}

for bucket in $(seq 0 $((TOTAL_BUCKETS - 1))); do
    echo "Starting training for bucket $bucket"
    bash submit_main.sh $SCHEDULE $MODE $bucket
done
