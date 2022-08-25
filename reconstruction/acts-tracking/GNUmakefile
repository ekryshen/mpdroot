TARGET = libMpdTpcTracker.so

TEST = test

SOURCES = ActsExamples/TGeoDetector/TGeoDetector.cxx \
          ActsExamples/TGeoDetector/TGeoITkModuleSplitter.cxx \
          MpdTpcConfig.cxx \
          MpdTpcDetector.cxx \
          MpdTpcDigitization.cxx \
          MpdTpcMagneticField.cxx \
          MpdTpcRunner.cxx \
          MpdTpcSpacePointMaking.cxx \
          MpdTpcTrackEstimation.cxx \
          MpdTpcTrackFinding.cxx \
          MpdTpcTrackSeeding.cxx

OBJECTS = $(SOURCES:%.cxx=%.o)

TEST_SOURCES = test.cxx
TEST_OBJECTS = $(TEST_SOURCES:%.cxx=%.o)

CC = g++

ACTS_HOME = $(DUBNA_HOME)/acts/git/build
FAIRSOFT_HOME = $(DUBNA_HOME)/install/FairSoft

INCLUDES = -I. -I$(FAIRSOFT_HOME)/include/root -I$(FAIRSOFT_HOME)/include

CFLAGS = -std=c++17 -Wall -Wextra -fPIC -g3
LFLAGS = -L$(ACTS_HOME)/lib \
         -L$(FAIRSOFT_HOME)/lib \
         -lActsCore -lActsPluginTGeo -lCore -lGeom -lboost_system -lboost_program_options

all: $(TARGET) $(TEST) install

$(TARGET): $(OBJECTS)
	$(CC) $^ -shared -o $@ $(LFLAGS)

$(TEST): $(TEST_OBJECTS) $(TARGET)
	$(CC) $^ -o $@ $(LFLAGS)

install: $(TARGET)
	cp $(TARGET) $(FAIRSOFT_HOME)/lib # :-)

%.o: %.cxx
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: clean run

run:
	./$(TEST)

debug:
	gdb --args ./$(TEST)

view:
	meshlab tpc_acts_tracking.obj

clean:
	rm -rf $(TARGET) $(TEST) $(OBJECTS) $(TEST_OBJECTS)
