# GNUmakefile 

PROG = reco

C++         = g++
DEBUG_OPTS  = -g -Wall
#ROOTCINT    = rootcint

FSDIR = /jinr/fairsoft/install
FRDIR = /jinr/fairroot/install
MRDIR = /home/fomenko/work/proj/mpd/mpdroot
#MIDIR = /home/fomenko/work/proj/mpd/tpc-study

ALL_LIB = $(shell root-config --glibs) -lMinuit -lXMLIO -lXMLParser
ALL_INC = $(shell root-config --cflags) -I${FSDIR}/include -I${FSDIR}/include/vmc -I${FRDIR}/include 

ALL_LIB += -L${FRDIR}/lib -L${FSDIR}/lib -L${MRDIR}/build/lib -lMpdBase -ltpc -lMpdField -lMpdMCStack -lMpdGen -lBase -lParBase -lFairLogger
ALL_INC += -I/usr/include/eigen3 -I${MRDIR}/mpdfield -I${MRDIR}/tpc -I${MRDIR}/mcstack -I${MRDIR}/generators -I${MRDIR}/kalman

all: ${PROG}

${PROG}: ${PROG}.C
	$(C++) $(DEBUG_OPTS) $(ALL_INC)  $^ -o $@ ${ALL_LIB}

#dict:
#	${ROOTCINT} -f ${PROG}_d.cxx -c ${PROG}.h LinkDefs.h

clean: 
	rm -rf ${PROG} ${PROG}_d.*

