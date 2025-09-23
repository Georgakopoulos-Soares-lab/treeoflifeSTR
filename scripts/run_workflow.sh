#!/bin/bash

SCHEDULE=$1
PATTERN=$2
TOTAL_BUCKETS=$3
CORES=$4


bash run_main.sh $SCHEDULE $PATTERN $TOTAL_BUCKETS
bash run_snake.sh $CORES
