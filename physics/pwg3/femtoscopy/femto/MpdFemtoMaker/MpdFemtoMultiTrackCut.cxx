//
// Allows to add several cuts to the list
//

// MpdFemtoMaker headers
// Infrastructure
#include "MpdFemtoMultiTrackCut.h"
#include "MpdFemtoTrackCutCollection.h"

ClassImp(MpdFemtoMultiTrackCut)

//_________________
MpdFemtoMultiTrackCut::MpdFemtoMultiTrackCut() {
  // Default constructor
  mCutCollection = new MpdFemtoTrackCutCollection;
}

//_________________
MpdFemtoMultiTrackCut::MpdFemtoMultiTrackCut(const MpdFemtoMultiTrackCut& copy) : MpdFemtoBaseTrackCut() {
  // Copy constructor
  MpdFemtoTrackCutIterator iter;
  MpdFemtoBaseTrackCut *t;
  for (iter = copy.mCutCollection->begin(); iter != copy.mCutCollection->end(); iter++) {
    t = (*iter)->clone();
    if (t) {
      mCutCollection->push_back(t);
    }
  } //for (iter=copy.mCutCollection->begin(); iter!=copy.mCutCollection->end(); iter++)
}

//_________________
MpdFemtoMultiTrackCut& MpdFemtoMultiTrackCut::operator=(const MpdFemtoMultiTrackCut& copy) {
  if (!mCutCollection) {
    MpdFemtoMultiTrackCut();
  }

  if (this != &copy) {
    MpdFemtoTrackCutIterator iter;
    MpdFemtoBaseTrackCut *t;
    for (iter = copy.mCutCollection->begin(); iter != copy.mCutCollection->end(); iter++) {
      t = (*iter)->clone();
      if (t) {
	mCutCollection->push_back(t);
      } //if (t)
    } // for (iter=copy.mCutCollection->begin(); iter!=copy.mCutCollection->end(); iter++)
  } //if ( this != &copy )

  return *this;
}

//_________________
MpdFemtoMultiTrackCut::~MpdFemtoMultiTrackCut() {
  if (mCutCollection) {
    MpdFemtoTrackCutIterator iter;
    for (iter = mCutCollection->begin(); iter != mCutCollection->end(); iter++) {
      delete *iter;
    }
  }
  delete mCutCollection;
}

//_________________
MpdFemtoString MpdFemtoMultiTrackCut::report() {
  MpdFemtoString temp;
  temp = "\n MpdFemtoMultiTrackCut::report() reporting ";
  temp += (int) mCutCollection->size();
  temp += " track cuts \n";
  MpdFemtoTrackCutIterator iter;
  for (iter = mCutCollection->begin(); iter != mCutCollection->end(); iter++) {
    temp += (*iter)->report();
  }
  return temp;
}

//_________________
bool MpdFemtoMultiTrackCut::pass(const MpdFemtoTrack* track) {
  bool temp = 0;
  MpdFemtoTrackCutIterator iter;
  for (iter = mCutCollection->begin(); iter != mCutCollection->end(); iter++) {
    temp = temp || (*iter)->pass(track);
  }
  return temp;
}

//_________________
void MpdFemtoMultiTrackCut::eventBegin(const MpdFemtoEvent* ev) {
  MpdFemtoTrackCutIterator iter;
  for (iter = mCutCollection->begin(); iter != mCutCollection->end(); iter++) {
    (*iter)->eventBegin(ev);
  }
}

//_________________
void MpdFemtoMultiTrackCut::eventEnd(const MpdFemtoEvent* ev) {
  MpdFemtoTrackCutIterator iter;
  for (iter = mCutCollection->begin(); iter != mCutCollection->end(); iter++) {
    (*iter)->eventEnd(ev);
  }
}

//_________________
void MpdFemtoMultiTrackCut::addTrackCut(MpdFemtoBaseTrackCut* cut) {
  mCutCollection->push_back(cut);
}

//_________________
MpdFemtoMultiTrackCut* MpdFemtoMultiTrackCut::clone() {
  MpdFemtoMultiTrackCut* cut = this->clone();
  return cut;
}
