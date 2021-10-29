#!/bin/bash

filelist_pile=/SMASH_txt/chunks/

for k in {0..199}
do
qsub -V setup.sh  $filelist_pile$k.lis
done


