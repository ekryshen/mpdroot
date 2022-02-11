#!/bin/bash

#
#$ -wd /scratch1/maldonado/Soft/Centrality/ana
#$ -cwd
#$ -N run_centrality_fit
#$ -q all.q
#$ -l h_rt=8:00:00
#$ -l s_rt=8:00:00
#
#$ -o /scratch1/maldonado/Soft/Centrality/ana/LOGS
#$ -e /scratch1/maldonado/Soft/Centrality/ana/LOGS
#

START_DIR=$PWD

cd /scratch1/maldonado/Soft/Centrality/ana

# Sourcing ROOT
source /cvmfs/nica.jinr.ru/sw/os/login.sh
module add mpddev
export MPDROOT=/lhep/users/maldonado/mySoft/mpd
source $MPDROOT/config/env.sh

cp ${START_DIR}/runanalysis$JOBID.C .

root -b -q runanalysis$JOBID.C
cd ..
