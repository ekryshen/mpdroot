#!/bin/bash
# export MPDROOT=$(MPDROOT-@CMAKE_INSTALL_PREFIX@) # if defined MPDROOT use MPDROOT, else use CMAKE_INSTALL_PREFIX
export MPDROOT="$(dirname $(dirname $(readlink -m ${BASH_ARGV[0]})))"
export LD_LIBRARY_PATH=$MPDROOT/lib:$LD_LIBRARY_PATH
export PATH=$MPDROOT:$PATH
export ROOT_INCLUDE_PATH=$MPDROOT/macros/common:$MPDROOT/include:$ROOT_INCLUDE_PATH
export GEOMPATH=$MPDROOT/geometry
export CONFIG_DIR=$MPDROOT/gconfig

export VMCWORKDIR=$MPDROOT # temporary patch until VMCWORKDIR will be removed
 # export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:"@MY_DYLD_LIBRARY_PATH@"
