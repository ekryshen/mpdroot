#!/bin/bash

#
# Assuming that this macro is placed in its original location
# (in the MpdMiniEvent directory) its path can be used to set
# the MpdMiniDst enviroment
#

filnam=$BASH_SOURCE            # filename of this (sourced) macro
absnam=`readlink -f $filnam`   # absolute filename
pather=`dirname $absnam`       # path to the config directory

export MINIDST=$pather

export PATH=${MINIDST}':'${PATH}
export LD_LIBRARY_PATH=${MINIDST}':'${LD_LIBRARY_PATH}

echo MINIDST = $MINIDST
