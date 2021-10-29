#!/bin/bash

export ENERGY=09GeV
export INDIR=dcap://lxse-dc01.jinr.ru:22126//pnfs/jinr.ru/data/mpd/data4mpd/data/
export MACDIR=/opt/exp_soft/mpd/geger/new_prod/$ENERGY/
export OUTDIR=dcap://lxse-dc01.jinr.ru:22126//pnfs/jinr.ru/data/mpd/data4mpd/dst/$ENERGY/
export MACFILE="reco.C"
export NEVENTS=200

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/exp_soft/d-cache/dcap/lib64
export LD_PRELOAD=libpdcap.so
export DCACHE_IO_TUNNEL=libgssTunnel.so

for INFILE in ` ls $INDIR$ENERGY ` ; do
  export INFILE
  export OUTFILE=dst_$INFILE
  qsub -v MACDIR,MACFILE,INDIR,INFILE,OUTDIR,OUTFILE,NEVENTS,ENERGY reco_from_dcache.sh
done

