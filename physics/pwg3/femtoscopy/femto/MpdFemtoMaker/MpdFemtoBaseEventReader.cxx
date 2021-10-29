//
// The pure virtual base class for femto event readers
//

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseEventCut.h"
#include "MpdFemtoBaseTrackCut.h"
#include "MpdFemtoBaseV0Cut.h"
#include "MpdFemtoBaseXiCut.h"
#include "MpdFemtoBaseKinkCut.h"
#include "MpdFemtoBaseEventReader.h"
// Infrastructure
#include "MpdFemtoEvent.h"

//_________________
MpdFemtoBaseEventReader::MpdFemtoBaseEventReader() :
mEventCut(nullptr),
mTrackCut(nullptr),
mV0Cut(nullptr),
mXiCut(nullptr),
mKinkCut(nullptr),
mReaderStatus(0),
mDebug(1) {
  /* empty */
}

//_________________
MpdFemtoBaseEventReader::MpdFemtoBaseEventReader(const MpdFemtoBaseEventReader &copy) :
mEventCut(copy.mEventCut),
mTrackCut(copy.mTrackCut),
mV0Cut(copy.mV0Cut),
mXiCut(copy.mXiCut),
mKinkCut(copy.mKinkCut),
mReaderStatus(copy.mReaderStatus),
mDebug(copy.mDebug) {
  /* empty */
}

//_________________
MpdFemtoBaseEventReader& MpdFemtoBaseEventReader::operator=(const MpdFemtoBaseEventReader& aReader) {
  if (this != &aReader) {
    mEventCut = aReader.mEventCut;
    mTrackCut = aReader.mTrackCut;
    mV0Cut = aReader.mV0Cut;
    mXiCut = aReader.mXiCut;
    mKinkCut = aReader.mKinkCut;
    mReaderStatus = aReader.mReaderStatus;
    mDebug = aReader.mDebug;
  }

  return *this;
}

//_________________
MpdFemtoString MpdFemtoBaseEventReader::report() {

  MpdFemtoString temp = "\n This is the base class MpdFemtoBaseEventReader reporting";
  temp += "\n---> EventCuts in Reader: ";
  if (mEventCut) {
    temp += mEventCut->report();
  } else {
    temp += "NONE";
  }
  temp += "\n---> TrackCuts in Reader: ";
  if (mTrackCut) {
    temp += mTrackCut->report();
  } else {
    temp += "NONE";
  }
  temp += "\n---> V0Cuts in Reader: ";
  if (mV0Cut) {
    temp += mV0Cut->report();
  } else {
    temp += "NONE";
  }
  temp += "\n---> XiCuts in Reader: ";
  if (mXiCut) {
    temp += mXiCut->report();
  } else {
    temp += "NONE";
  }
  temp += "\n---> KinkCuts in Reader: ";
  if (mKinkCut) {
    temp += mKinkCut->report();
  } else {
    temp += "NONE";
  }
  temp += "\n";
  return temp;
}

//_________________
int MpdFemtoBaseEventReader::init(const char* /* ReadWrite */, MpdFemtoString& /* Message */) {
  std::cout << "do-nothing MpdFemtoBaseEventReader::Init()" << std::endl;
  return 0;
}

ClassImp(MpdFemtoBaseEventReader)
