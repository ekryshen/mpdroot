#!/bin/bash
#
# qsub -p +1023 pbsrun.sh
# qsub -l walltime=5:30,cput=5:30,mem=256mb pbsrun.sh
# qstat | grep $USER
# qstat -f 997081.lxbsrv01.jinr.ru
#
# qdel 997081.lxbsrv01.jinr.ru
#
export INFILE

export DATDIR=$INDIR$ENERGY/
export DATFILE=$INFILE

export OUTDIR=$MACDIR
export OUTFILE=a$INFILE

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/exp_soft/d-cache/dcap/lib64
export LD_PRELOAD=libpdcap.so
export DCACHE_IO_TUNNEL=libgssTunnel.so

cd $TMPDIR
dccp $DATDIR$DATFILE .
cp  $MACDIR$MACFILE .

. /opt/exp_soft/mpd/fairsoft_4.1.2/20100115/SetEnv_4.1.2.sh

export LD_HOME=/opt/exp_soft/mpd/mpdroot
export HOME=/opt/exp_soft/mpd/trunk

export VMCWORKDIR=$HOME
export LD_LIBRARY_PATH=$LD_HOME/lib:$LD_LIBRARY_PATH

root -q -b $MACFILE\(\"$INFILE\"\)

cp $OUTFILE $OUTDIR
rm $INFILE
