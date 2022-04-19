# Source files (*.C)
#SRCS=runMC.C
SRCS=reco.C

# Compiler
CC=g++

OBJECTS=$(patsubst %.C, %, $(SRCS))

ROOTCONFIG := root-config

CFLAGS := $(shell $(ROOTCONFIG) --cflags)
CFLAGS += -I$(VMCWORKDIR)/include
CFLAGS += -I${FAIRROOTPATH}/include
CFLAGS += -I$(SIMPATH)/include
CFLAGS += -I$(ROOT_INCLUDE_PATH)
CFLAGS += -I$(ROOT_INCLUDE_DIR)
ROOT_PATH = -I$(ROOT_INCLUDE_PATH) -I$(ROOT_INCLUDE_DIR)
CFLAGS += -g -O0

LDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
LDFLAGS += -lEG $(shell $(ROOTCONFIG) --glibs)
LDFLAGS += -L$(VMCWORKDIR)/build/lib
LDFLAGS += -L$(ROOT_LIBRARY_DIR) -L$(FAIRROOTPATH)/lib
LDFLAGS += -lKalman -lLHETrack -lMpdPid -lMpdDst
LDFLAGS += -ltpc -lTof -lEtof -lEmc -lZdc -lSts -lFfd
LDFLAGS += -lMpdBase -lPassive -lMpdField -lMpdGen
LDFLAGS += -lBase -lGeoBase -lParBase -lTrkBase -lGeane -lFairTools


all:
	$(OBJECTS)

$(OBJECTS): % : $(SRCS)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -vf $(basename $(SRCS)) *.o

