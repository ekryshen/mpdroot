#!/bin/bash

# Governor script for running build tests in Gitlab's CI pipeline
# Written for MPDRoot @JINR Dubna
# First commit Slavomir Hnatic, 12.2021

# NOTE: This governor is not easy to modify
#       Do not complain as other frameworks are more difficult
#       (if in doubt have a look at Boost::test)

##########################################################################################################
########################################### Functionality ################################################
###                                                                                                    ###
###  1. INPUT PARAMETERS                                                                               ###
###                                                                                                    ###
###     -t ...optional. If this option is passed then old type build tests are executed. Otherwise     ###
###           the script defaults to running build tests for the new build system.                     ###
###     -r ...REPORT_NAME. Assigns REPORT_FILENAME=$REPORT_NAME.xml (default = report.xml)             ###
###                                                                                                    ###
###  2. OUTPUT                                                                                         ###
###                                                                                                    ###
###     written into $REPORT_FILENAME (json xml format), which is then picked by Gitlab to generate    ###
###     graphical output in "Tests" section of the Gitlab's pipeline                                   ###
###                                                                                                    ###
##########################################################################################################
######################################### Adding new tests ###############################################
###                                                                                                    ###
###  1.  Add new pair in function assign_vmc_generator_pairs, increment NO_VMC_GEN_PAIRS               ###
###                                                                                                    ###
###  1.1 Add new test into vmc_generator test suite (inside "for" cycle in function generate_tests)    ###
###                                                                                                    ###
###  2.  Add new test suite (after the "for" cycle in function generate_tests)                         ###
###       - you might need to create workdir for your newly created tests (in function run_tests)      ###
###                                                                                                    ###
###      ALWAYS MAKE SURE YOU KNOW WHAT YOU ARE DOING, IF IN DOUBT ASK SUPPORT TO ADD TEST FOR YOU     ###
###      BY CREATING A GITLAB ISSUE ON https://git.jinr.ru/nica/mpdroot/-/issues)                      ###
###                                                                                                    ###
##########################################################################################################



main() {

 # assign parameters
 MPDROOT_MACROS=$CI_PROJECT_DIR/macros
 while getopts t:r: flag
 do
    case "${flag}" in
         t) OLD_BUILD_TYPE=${OPTARG};;
         r) REPORT_NAME=${OPTARG};;
    esac
 done
 [[ $REPORT_NAME ]] && REPORT_FILENAME="$REPORT_NAME.xml" || REPORT_FILENAME="report.xml"

 if [[ $OLD_BUILD_TYPE ]]; then
  [[ $OS_TYPE == "centos7" ]] && source /opt/rh/devtoolset-7/enable
  source build/config.sh
  generate_old_tests
 else
  assign_vmc_generator_pairs
  generate_tests
 fi

 run_tests
 write_test_report

}


function assign_vmc_generator_pairs() {

  # definition of pairs (VMC, GENERATOR) for which tests are generated.
  # MPD generators: BOX FLUID HSD ION LAQGSM MCDST PART SMASH UNIGEN URQMD VHLLE
  # MPD vmcs : GEANT3 GEANT4
  # Each pair is one test suite for which vmc,generator test template is executed
  NO_VMC_GEN_PAIRS=4
  INFILE_DEFAULT="auau.09gev.mbias.98k.ftn14"
  vmc[0]="GEANT4"; generator[0]="BOX"; inFile[0]=$INFILE_DEFAULT;
  vmc[1]="GEANT4"; generator[1]="UNIGEN"; inFile[1]="$CI_PROJECT_DIR/input/tests/dcmqgsm_bibi_9.2gev_local_1.mcini.root";
  vmc[2]="GEANT3"; generator[2]="BOX";    inFile[2]=$INFILE_DEFAULT;
  vmc[3]="GEANT3"; generator[3]="UNIGEN"; inFile[3]="$CI_PROJECT_DIR/input/tests/dcmqgsm_bibi_9.2gev_local_1.mcini.root";
}


function generate_tests() {

 IFS=
 echo "Generating default (new) style tests"

 #iterate over vmc/gen pairs
 TEST_COUNTER=0
 TESTS_IN_SUITE=3
 for i in "${!vmc[@]}"; do
    tests[$TEST_COUNTER]='root -b -q -l '"'"'$MPDROOT_MACROS/common/runMC.C(EGenerators::'"${generator[$i]}"',EVMCType::'"${vmc[$i]}"',1,"'"${inFile[$i]}"'","evetest.root",0,1)'"'"' | tee output_'"$TEST_COUNTER"'.txt'
    names[$TEST_COUNTER]="runMC"
    suitename[$TEST_COUNTER]="${vmc[$i]} ${generator[$i]}"
    ((TEST_COUNTER++))
    tests[$TEST_COUNTER]='root -b -q -l '"'"'$MPDROOT_MACROS/common/runReco.C()'"'"' | tee output_'"$TEST_COUNTER"'.txt'
    names[$TEST_COUNTER]="reco"
    suitename[$TEST_COUNTER]="${vmc[$i]} ${generator[$i]}"
    ((TEST_COUNTER++))
    tests[$TEST_COUNTER]='root -b -q -l '"'"'$MPDROOT_MACROS/common/readDST.C("mpddst.root")'"'"' | tee output_'"$TEST_COUNTER"'.txt'
    names[$TEST_COUNTER]="readDST"
    suitename[$TEST_COUNTER]="${vmc[$i]} ${generator[$i]}"
    ((TEST_COUNTER++))
 done
 TESTS_VMC_GEN=$TEST_COUNTER

}


function generate_old_tests() {

  IFS=
  echo "Generating old style tests"

  tests[0]='root -b -q -l '"'"'$CI_PROJECT_DIR/macro/mpd/runMC.C("auau.09gev.mbias.98k.ftn14","$VMCWORKDIR/macro/mpd/evetest.root",0,2)'"'"' | tee output_0.txt'
  names[0]="runMC"
  suitename[0]="Old build"
  tests[1]='root -b -q -l '"'"'$CI_PROJECT_DIR/macro/mpd/reco.C'"'"' | tee output_1.txt'
  names[1]="reco"
  suitename[1]=${suitename[0]}

}


function run_tests() {

 echo "Running tests"

 INDEX=0
 for test in ${tests[@]}; do

  if [[ ! $OLD_BUILD_TYPE && $(( $INDEX % $TESTS_IN_SUITE )) == 0 && (( $INDEX < $TESTS_VMC_GEN )) ]]; then
    # enter workdir
   INDEX_VMC=$(( $INDEX / $TESTS_IN_SUITE ))
   cd $CI_PROJECT_DIR && mkdir -p ${vmc[$INDEX_VMC]}/${generator[$INDEX_VMC]} && cd ${vmc[$INDEX_VMC]}/${generator[$INDEX_VMC]}/
  fi

  echo $test && ls -l
  ( time -p eval $test ) 2> err.txt
  timings[INDEX]=$(cat err.txt | grep real | sed "s/^[^ ]* //" )
  echo "Test time is ${timings[$INDEX]} seconds."

  if grep -q "Macro finished successfully." output_$INDEX.txt; then
   echo "${names[$INDEX]} test passed."
   result[$INDEX]=pass
   err_output[$INDEX]=""
  else
   echo "${names[$INDEX]} test failed"
   result[$INDEX]=fail
   err_output[$INDEX]=$(tail -n 25 output_$INDEX.txt)
  fi

  ((INDEX++))
 done

}


function write_test_report() {

 echo "Writing report into file"
 cd $CI_PROJECT_DIR

 echo "<testsuite name=\"Build test suite\" tests=\"$INDEX\">" > $REPORT_FILENAME

 INDEX=0
 for test in ${tests[@]}; do
  echo "    <testcase classname=\"${suitename[$INDEX]}\" name=\"${names[$INDEX]}\" time=\"${timings[$INDEX]}\">" >> $REPORT_FILENAME
  if [[ ${result[$INDEX]} == "fail" ]]; then
   echo "        <failure type=\"MacroDidNotFinish\"> ${err_output[$INDEX]}  </failure>" >> $REPORT_FILENAME
  fi
  echo "    </testcase>" >> $REPORT_FILENAME
  ((INDEX++))
 done

 echo "</testsuite>" >> $REPORT_FILENAME

}


main "$@"

