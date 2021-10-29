#!/bin/bash

export ENERGY=09GeV
export INDIR=dcap://lxse-dc01.jinr.ru:22126//pnfs/jinr.ru/data/mpd/data4mpd/dst/
export MACDIR=/opt/exp_soft/mpd/geger/new_prod/$ENERGY/
export MACFILE="anaDST.C"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/exp_soft/d-cache/dcap/lib64
export LD_PRELOAD=libpdcap.so
export DCACHE_IO_TUNNEL=libgssTunnel.so

for INFILE in ` ls $INDIR$ENERGY ` ; do
  export INFILE,
  qsub -v INDIR,INFILE,ENERGY,MACDIR,MACFILE ana_from_dcache.sh
done

