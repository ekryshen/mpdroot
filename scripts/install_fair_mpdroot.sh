#! /bin/bash

# Script to install FairSoft vApr21, FairRoot v18.6.4 & Mpdroot
# Written for MPDRoot @JINR Dubna
# First commit Slavomir Hnatic, 12.2021

#################################################################################
############################## Functionality ####################################
###                                                                           ###
###  1.INPUT PARAMETERS (default values in scopes)                            ###
###                                                                           ###
###    -j, --cpu_threads (8)         ...no of threads used in compilation     ###
###    -u, --user_type (regular)     ...developer or regular user             ###
###    -f, --fair_dir (~/fairSuite)  ...FairSoft/FairRoot                     ###
###                                     installation directory                ###
###    -m, --mpdroot_dir (~/mpd)     ...MpdRoot installation directory        ###
###    -s, --skip_fair               ...if set to yes                         ###
###                                     only mpdroot is installed             ###
###    -w, --fair_workdir (fair_build)  ...work directory name to build       ###
###                                        Fair Suite                         ###
###    -r, --mpdroot_repodir (mpdroot)  ...directory name where mpdroot       ###
###                                        repository will be cloned          ###
###                                                                           ###
###                                                                           ###
###  2.EXAMPLES (running the script from your $HOMEDIR)                       ###
###                                                                           ###
###    ./install_fair_mpdroot.sh -j 16 -u developer -s yes                    ###
###    installs mpdroot if Fair suite is already installed in default dir     ###
###    for the developer user with 16 threads used in compilation             ###
###                                                                           ###
###    ./install_fair_mpdroot.sh -m ~/mpdDec21 -r ~/mpdrootDec21              ###
###    install Fair Suite with default parameters and install mpdroot into    ###
###    ~/mpdDec21 directory and clone repository into ~/mpdrootDec21 dir      ###
###                                                                           ###
###    ./install_fair_mpdroot.sh                                              ###
###    installs Fair suite and mpdroot with default parameters                ###
###                                                                           ###
###  3.All sources/build of Fair suite are in created $fair_workdir dir       ###
###    mpdroot build is located in $mpdroot_repodir/build                     ###
###                                                                           ###
#################################################################################


main() {

 # exit when any command fails
 set -e

 # init/read parameter values
 cpu_threads="8"
 user_type="regular"
 fair_dir=~/fairSuite
 mpdroot_dir=~/mpd
 skip_fair=""
 fair_workdir=fair_build
 mpdroot_repodir=mpdroot

 while getopts j:u:f:m:s:w:r: flag
 do
    case "${flag}" in
         j) cpu_threads=${OPTARG};;
         u) user_type=${OPTARG};;
         f) fair_dir=${OPTARG};;
         m) mpdroot_dir=${OPTARG};;
         s) skip_fair=${OPTARG};;
         w) fair_workdir=${OPTARG};;
         r) mpdroot_repodir=${OPTARG};;
    esac
 done

 init_dir=$(pwd)

 # install Fair Suite if needed
 if [[ "$skip_fair" != "yes" ]]; then
  install_fairsuite
 fi

 # install mpdroot
 install_mpdroot

 # write out the installation summary
 generate_summary

}


function install_fairsuite() {

 # FairSoft build and install
 mkdir $fair_workdir && cd $fair_workdir
 git clone -b apr21p1 https://github.com/FairRootGroup/FairSoft
 cd FairSoft && mkdir build && cd build
 cmake -DCMAKE_INSTALL_PREFIX=$fair_dir/FairSoft -DGEANT4MT=OFF ..
 make -j $cpu_threads
 sed -i '/#include "UserDefaults.h"/a #include <thread>' Source/dds/dds-info/src/main.cpp
 make -j $cpu_threads
 sed -i '/#include <vector>/a #include <thread>' Source/fairmq/fairmq/sdk/DDSSession.cxx
 sed -i 's/  int dummy;/  int dummy = 0;/g' Source/fairmq/extern/googletest/googletest/src/gtest-death-test.cc
 sed -i 's/  bool result;/  bool result = false;/g' Source/fairmq/extern/googletest/googletest/src/gtest-death-test.cc
 make -j $cpu_threads
 cd $fair_dir/FairSoft/lib
 ln -s libpythia6.so libPythia6.so
 export SIMPATH=$fair_dir/FairSoft

 # FairRoot build and install
 cd $init_dir/$fair_workdir
 git clone -b v18.6.4 https://github.com/FairRootGroup/FairRoot.git
 cd FairRoot && mkdir build && cd build
 cmake -DCMAKE_INSTALL_PREFIX=$fair_dir/FairRoot ..
 make -j $cpu_threads
 make install

}

function install_mpdroot() {

 cd $init_dir
 if [ "$user_type" == "regular" ]; then
  git clone -b dev --recursive https://git.jinr.ru/nica/mpdroot.git $mpdroot_repodir
 else
  if [ "$user_type" == "developer" ]; then
   git clone -b dev --recursive git@git.jinr.ru:nica/mpdroot.git $mpdroot_repodir
   # add clang pre-commit hook
   cp $mpdroot_repodir/scripts/pre-commit $mpdroot_repodir/.git/hooks/
   chmod u+x $mpdroot_repodir/.git/hooks/pre-commit
  else
   echo -e "Only regular or developer user types allowed! Terminating...\n"
   exit 1
  fi
 fi

 # linking with local Fair suite installation
 cd $mpdroot_repodir
 sed -i "2i fair_dir=$fair_dir" SetEnv.sh
 sed -i 's/export SIMPATH=\/opt\/fairsoft\/install/export SIMPATH=$fair_dir\/FairSoft/g' SetEnv.sh
 sed -i 's/export FAIRROOTPATH=\/opt\/fairroot\/install/export FAIRROOTPATH=$fair_dir\/FairRoot/g' SetEnv.sh
 source SetEnv.sh

 # prevent SetEnv.sh from being overwritten by future git synchronizations
 git rm --cached SetEnv.sh

 # create & enter local build directory. run cmake config
 mkdir build && cd build && cmake -DCMAKE_INSTALL_PREFIX=$mpdroot_dir ..

 # build
 make -j $cpu_threads

 # source new configuration
 source config.sh

 # install
 make install

}


function generate_summary() {

echo -e "\n   INSTALLATION FINISHED SUCCESSFULLY
\n
 Fair Suite is installed in $fair_dir
 MpdRoot is installed in $mpdroot_dir
 Fair Suite sources and build are located in $init_dir/$fair_workdir
 MpdRoot build is located in $init_dir/$mpdroot_repodir/build
 MpdRoot cloned repository is located in $init_dir/$mpdroot_repodir
\n
 When using MPDRoot you must in each terminal load the environment by:
 source $mpdroot_repodir/SetEnv.sh && source $mpdroot_repodir/build/config.sh
\n
 You can add this line at the end of your ~/.bashrc file to do it automatically
\n
"

}


main "$@"

