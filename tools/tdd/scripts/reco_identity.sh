#!/bin/bash

# TDD script to test reconstruction for identical results
# Written for MPDRoot @JINR Dubna
# First commit Slavomir Hnatic, 5.2022

##########################################################################################################
########################################### Functionality ################################################
###                                                                                                    ###
###  1. INTRO                                                                                          ###
###     ... run in a separate directory from your local repo clone (a friendly suggestion)             ###
###     ... $INFILE, $OUTFILE, $OUTFILE_DIR are then stored in this directory                          ###
###                                                                                                    ###
###  2. INPUT PARAMETERS                                                                               ###
###                                                                                                    ###
###     any parameter ... generate $INFILE and reference output file $OUTFILE_REF                      ###
###     no parameter  ... generate $OUTFILE from existing $INFILE and compare with                     ###
###                       $OUTFILE_REF for identity                                                    ###
###                                                                                                    ###
##########################################################################################################


main() {

 # exit when any command fails
 set -e

 # check and set basic variables
 if [[ -z $MPDROOT ]]; then
     echo "Please set the MPDROOT variable (installation directory) !"
     exit 1
 fi

 RUNMC=$MPDROOT/macros/common/runMC.C
 RUNRECO=$MPDROOT/macros/common/runReco.C

 INFILE=evetest.root

 TMPFILE=tmp.txt
 OUTFILE=mpddst.txt
 OUTFILE_REF=mpddst_ref.txt

 # generate output with comparison to reference || generate reference output
 [[ -z $1 ]] &&  generate_comparison || generate_reference

 file_cleanup

}


function generate_comparison() {
 if [[ ! $INFILE || ! $OUTFILE_REF ]]; then
  echo -e "Reference data do not exist.\nRun the script with generate reference data option first."
  exit 1
 fi
 root -l -b -q "$RUNRECO" > $TMPFILE
 filter > $OUTFILE

 echo ""
 echo -e "\033[32mIDENTITY CHECK RESULT:\033[0m"
 diff -s $OUTFILE $OUTFILE_REF
 echo ""
}

function generate_reference() {
 touch $INFILE && rm $INFILE #remove $INFILE if exists
 root -l -b -q "$RUNMC"
 root -l -b -q "$RUNRECO" > $TMPFILE
 filter > $OUTFILE_REF
}

function filter() {
 sed '1,/\[INFO\] FairRunAna/d' $TMPFILE | sed '/Digitizer work time/,$d' | sed '/User CPU time/d'
}

function file_cleanup() {
 THISFILE=`basename "$0"`
 find . -maxdepth 1 -type f \! \( -name "$THISFILE" -o -name "$INFILE" -o -name "$OUTFILE_REF" \) -delete
}

main "$@"


