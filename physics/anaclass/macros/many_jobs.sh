#!/bin/bash

for ((INDEX=340; INDEX < 345; INDEX++))
do
   sed -e "s/runanalysis/runanalysis$INDEX/; s/urqmd-BiBi-09.2GeV-mb-eos0-500-347/urqmd-BiBi-09.2GeV-mb-eos0-500-$INDEX/; s/outputana0.root/outputana$INDEX.root/" runanalysis.C > runanalysis$INDEX.C
qsub -v JOBID=$INDEX analysis.sh
done
