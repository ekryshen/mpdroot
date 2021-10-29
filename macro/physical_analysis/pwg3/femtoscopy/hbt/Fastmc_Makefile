CXX           =
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

CXX           = g++
CXXFLAGS      = -Wall -fPIC
LD            = g++
LDFLAGS       =
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS) -lEG
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

FASTMCO        = Fastmc_CF3d.o
UHKMS        = Fastmc_CF3d.C
 

              
TARGET	    = Fastmc_CF3d
#------------------------------------------------------------------------------
$(TARGET):       $(FASTMCO)
		$(LD)  -O2 $(LDFLAGS) $^ $(OutPutOpt) $@ $(LIBS)
		@echo "$@ done"
clean:
		@rm -f $(FASTMCO) $(TARGET)

%.o : %.cxx
	$(CXX) -O2 -I ../inc $(CXXFLAGS) -c $<
