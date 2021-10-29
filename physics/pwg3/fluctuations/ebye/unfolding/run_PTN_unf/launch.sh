#!/bin/bash

for n_subsamples in 1 5 10 20 50
do
for n_unf_iterations in 1 5 10 20 50
do
for n_bins in 5 10 15 20 50 75 90
do

echo ${n_subsamples}
echo ${n_unf_iterations}
echo ${n_bins}

qsub -V /run_PTN_unf/my_start.sh ${n_subsamples} ${n_unf_iterations} ${n_bins}

done
done
done
