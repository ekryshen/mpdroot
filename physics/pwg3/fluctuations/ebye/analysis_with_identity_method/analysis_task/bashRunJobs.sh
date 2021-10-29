#!/bin/bash

echo "#####  starting..."

step=20

for i in {0..9}; do
    startId=$(( $i*$step ))
    stopId=$(( ($i+1)*$step -1 ))
    echo $startId $stopId
    qsub sge_IA_job.qsub $startId $stopId
done
