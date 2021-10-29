#!/bin/bash
export ENERGY=04GeV
export INDIR=/opt/exp_soft/mpd/data4mpd/UrQMD/2.3/
if [ $ENERGY -eq 04GeV ]; then
    export INFILE=auau.04gev.0_3fm.10k.f14
elif [ $ENERGY -eq 07GeV ]; then
    export INFILE=auau.07gev.0-3fm.26.8k.f14
elif [ $ENERGY -eq 09GeV ]; then
    export INFILE=auau.09gev.0_3fm.41.5k.f14
elif [ $ENERGY -eq 11GeV ]; then
    export INFILE=auau.11gev.0-3fm.11k.f14
fi
if [ ! -f $INDIR$INFILE ]; then
  export INDIRTMP=dcap://lxse-dc01.jinr.ru:22126//pnfs/jinr.ru/data/mpd/data4mpd/UrQMD/2.3/
  dccp $INDIRTMP$INFILE $INDIR
fi

export INFILE=$INDIR$INFILE
export OUTDIR=dcap://lxse-dc01.jinr.ru:22126//pnfs/jinr.ru/data/mpd/data4mpd/data/$ENERGY/
export MACDIR=/opt/exp_soft/mpd/geger/new_prod/$ENERGY/
export MACFILE="runMC.C"

export OFILE=auau_04gev_0_3fm
export NFILES=50
export NEVENTS=200

i=0

while [ $i -lt $NFILES ]  
  do
  OUTFILE=$OFILE\_$i.root
  let "NSKIP = i * NEVENTS"
  export INFILE
  export OUTFILE
  export NSKIP
  qsub -v MACDIR,MACFILE,INDIR,INFILE,OUTDIR,OUTFILE,NEVENTS,NSKIP mc_dcache.sh
  let  "i = i + 1"
done

