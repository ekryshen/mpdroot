#!/bin/bash

export INDIR=dcap://lxse-dc01.jinr.ru:22126//pnfs/jinr.ru/data/mpd/data4mpd/dst/09GeV

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/exp_soft/d-cache/dcap/lib64
export LD_PRELOAD=libpdcap.so
export DCACHE_IO_TUNNEL=libgssTunnel.so

for INFILE in ` $INDIR ` ; do
  rm -f $INDIR$INFILE
done

