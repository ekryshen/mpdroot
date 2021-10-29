#!/bin/bash
#
# qsub -p +1023 pbsrun.sh
# qsub -l walltime=5:30,cput=5:30,mem=256mb pbsrun.sh
# qstat | grep ram
# qstat -f 997081.lxbsrv01.jinr.ru
#
# qdel 997081.lxbsrv01.jinr.ru
#

cd $MACDIR
#cd $TMPDIR
#cp $INDIR$INFILE $MACDIR$MACFILE .

. /opt/exp_soft/mpd/fairsoft_4.1.2/20100115/SetEnv_4.1.2.sh

export LD_HOME=/opt/exp_soft/mpd/geger/mpdroot
export HOME=/opt/exp_soft/mpd/geger/trunk

export VMCWORKDIR=$HOME
export LD_LIBRARY_PATH=/opt/exp_soft/d-cache/dcap/lib64:$LD_HOME/lib:$LD_LIBRARY_PATH
export LD_PRELOAD=libpdcap.so
export DCACHE_IO_TUNNEL=libgssTunnel.so

export INFILE
export OUTFILE
export NSKIP
export NEVENTS

root -q -b $MACFILE\(\"$INFILE\",\"$OUTFILE\",$NEVENTS,$NSKIP\)

dccp $OUTFILE $OUTDIR
rm -f $OUTFILE
