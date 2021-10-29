#!/bin/bash

n_subsamples=$1
n_unf_iterations=$2

export MAIN_DIR=/run_N_unf

export START_DIR=${PWD}
export OUT_DIR=${MAIN_DIR}/OUT
export TMP_DIR=${MAIN_DIR}/TMP
export OUT=${OUT_DIR}
export OUT_LOG=${OUT}/log
export OUT_FILE=${OUT}/files
export TMP=${TMP_DIR}/TMP_${JOB_ID}_${SGE_TASK_ID}
export LOG=${OUT_LOG}/JOB_${JOB_ID}_${SGE_TASK_ID}.log
touch $LOG

mkdir -p $OUT_LOG
mkdir -p $OUT_FILE
mkdir -p $TMP/${JOB_ID}
mkdir -p ${OUT_FILE}/${JOB_ID}
export out_location=${OUT_FILE}/${JOB_ID}
echo $out_location &>> $LOG
echo $impact &>> $LOG

. /*add_your_path*/mpdroot/SetEnv.sh &>> $LOG
. /*add_your_path*/mpdroot/build/config.sh &>> $LOG

cd $TMP/${JOB_ID}
echo $PWD &>> $LOG

cp /run_N_unf/runme.C .
cp /run_N_unf/final_unfolding_N.cc .

ls &>> $LOG

root -l -b -q 'runme.C('${n_subsamples}','${n_unf_iterations}')' &>> $LOG

echo "Job is done!" &>> $LOG



