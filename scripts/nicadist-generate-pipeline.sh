#!/bin/bash

# Script to generate MPDroot CI pipeline file for nicadist docker images
# Written for MPDRoot @JINR Dubna
# Authors: Jan Busa
# First commit: Jan Busa 12.2021

# Rationale: script can be used to generate the pipeline file
# for any combination of generators, virtual MCs and OSs

outFile=".gitlab-ci-nicadist.yml" # output name, do not change
targets=('c7' 'c8') # known targets
generators=('HADGEN') # 'BOX' 'FLUID' 'HADGEN' 'HSD' 'ION' 'LAQGSM' 'MCDST' 'PART' 'SMASH' 'UNIGEN' 'URQMD' 'VHLLE' # mpd generators
vmcs=('GEANT3') # 'GEANT3' 'GEANT4'  available virtual Monte Carlo

> $outFile # clear file, if present
for target in ${targets[@]}; do
cat >> $outFile << EOL
nicadist-${target}-build:
  image: registry.gitlab.com/ndmspc/user/${target}:alibuild-dev
  stage: build
  before_script:
    - source /etc/bashrc
    - export NPROC=\$(cat /proc/cpuinfo | grep processor | wc -l)
  script:
    - export MPDROOT=\$(readlink -m ${target}/mpdroot)
    - module add mpddev/latest
    - mkdir build
    - cd build
    - cmake ..
    - make install -j \$NPROC
  artifacts:
    paths:
      - ${target}/mpdroot/
    expire_in: 2h
  tags:
    - cvmfs
  only:
    - merge_requests
  allow_failure: true

EOL

  for vmc in ${vmcs[@]}; do
    for generator in ${generators[@]}; do
cat >> $outFile << EOL
nicadist-${target}-physics-${vmc}-${generator}:
  image: registry.gitlab.com/ndmspc/user/${target}:alibuild-dev
  stage: test_runMC
  dependencies:
    - nicadist-${target}-build
  before_script:
    - source /etc/bashrc
  script:
    - module add mpddev/latest
    - export MPDROOT=\$(readlink -m ${target}/mpdroot)
    - source \$MPDROOT/etc/env.sh
    - mkdir -p ${target}/${vmc}/${generator}
    - rm -rf ${target}/${vmc}/${generator}/*
    - cd ${target}/${vmc}/${generator}/
    - root -b -q -l '\$MPDROOT_MACROS/common/runMC.C(EGenerators::${generator},EVMCType::${vmc})'
    - root -b -q -l '\$MPDROOT_MACROS/common/runReco.C()'
    - root -b -q -l '\$MPDROOT_MACROS/common/readDST.C("mpddst.root")'
  retry: 1      # retry script if it failed (sometimes if fails due to some unexpected particle). Maximum 2 retries.
  tags:
    - cvmfs
  only:
    - merge_requests
  allow_failure: true

EOL
    done
  done
#  echo "" >> $outFile
done
