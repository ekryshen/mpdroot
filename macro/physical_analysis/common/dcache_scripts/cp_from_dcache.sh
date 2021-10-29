#!/bin/bash

export INDIR=/opt/exp_soft/mpd/geger/new_prod/09GeV/
export OUTDIR=dcap://lxse-dc01.jinr.ru:22126//pnfs/jinr.ru/data/mpd/data4mpd/dst/09GeV/

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/exp_soft/d-cache/dcap/lib64
export LD_PRELOAD=libpdcap.so
export DCACHE_IO_TUNNEL=libgssTunnel.so

for INFILE in ` ls $OUTDIR ` ; do
dccp $OUTDIR$INFILE $INDIR
done
