#!/bin/bash

for n_subsamples in 1 5 10
do
for n_unf_iterations in 1 5 10 20 50
do

echo ${n_subsamples}
echo ${n_unf_iterations}

qsub -V /run_NFNB_unf/my_start.sh ${n_subsamples} ${n_unf_iterations}

done
done
