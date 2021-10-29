//
// A handler for cut monitors
//

// C++ headers
#include <iterator>

// MpdFemtoMaker headers
#include "MpdFemtoCutMonitorHandler.h"

ClassImp(MpdFemtoCutMonitorHandler)

//_________________
MpdFemtoCutMonitorHandler::MpdFemtoCutMonitorHandler() :
mCollectionsEmpty(true),
mPassColl(nullptr),
mFailColl(nullptr) {

  // Default constructor
  mPassColl = new MpdFemtoCutMonitorCollection();
  mFailColl = new MpdFemtoCutMonitorCollection();
}

//_________________
MpdFemtoCutMonitorHandler::MpdFemtoCutMonitorHandler(const MpdFemtoCutMonitorHandler& copy) :
mCollectionsEmpty(copy.mCollectionsEmpty),
mPassColl(nullptr),
mFailColl(nullptr) {

  // Copy constructor
  mPassColl = new MpdFemtoCutMonitorCollection(copy.mPassColl->begin(), copy.mPassColl->end());
  mFailColl = new MpdFemtoCutMonitorCollection(copy.mFailColl->begin(), copy.mFailColl->end());
}

//_________________
MpdFemtoCutMonitorHandler& MpdFemtoCutMonitorHandler::operator=(const MpdFemtoCutMonitorHandler& copy) {

  // Assignment operator
  if (this != &copy) {

    if (mPassColl) {
      mPassColl->clear();
      mPassColl->insert(mPassColl->begin(), copy.mPassColl->begin(), copy.mPassColl->end());
    } else {
      mPassColl = new MpdFemtoCutMonitorCollection(copy.mPassColl->begin(), copy.mPassColl->end());
    }

    if (mFailColl) {
      mFailColl->clear();
      mFailColl->insert(mFailColl->begin(), copy.mFailColl->begin(), copy.mFailColl->end());
    } else {
      mFailColl = new MpdFemtoCutMonitorCollection(copy.mFailColl->begin(), copy.mFailColl->end());
    }
  } // if ( this != &copy )

  return *this;
}

//_________________
MpdFemtoCutMonitorHandler::~MpdFemtoCutMonitorHandler() {

  // Default destructor
  delete mPassColl;
  delete mFailColl;
  // if (mPassColl) { delete mPassColl; mPassColl = nullptr; }
  // if (mFailColl) { delete mFailColl; mFailColl = nullptr; }
}

//_________________
void MpdFemtoCutMonitorHandler::fillCutMonitor(const MpdFemtoEvent* event, bool pass) {

  // Fill event cut monitors
  if (mCollectionsEmpty) return;

  MpdFemtoCutMonitorIterator iter;
  MpdFemtoBaseCutMonitor* CM;

  if (pass) {
    for (iter = mPassColl->begin(); iter != mPassColl->end(); iter++) {
      CM = *iter;
      CM->fill(event);
    }
  }//if ( pass)
  else {
    for (iter = mFailColl->begin(); iter != mFailColl->end(); iter++) {
      CM = *iter;
      CM->fill(event);
    }
  } //else
}

//_________________
void MpdFemtoCutMonitorHandler::fillCutMonitor(const MpdFemtoTrack* track, bool pass) {

  // Fill track cut monitors
  if (mCollectionsEmpty) return;

  MpdFemtoCutMonitorIterator iter;
  MpdFemtoBaseCutMonitor* CM;

  if (pass) {
    for (iter = mPassColl->begin(); iter != mPassColl->end(); iter++) {
      CM = *iter;
      CM->fill(track);
    }
  } else {
    for (iter = mFailColl->begin(); iter != mFailColl->end(); iter++) {
      CM = *iter;
      CM->fill(track);
    }
  }
}

//_________________
void MpdFemtoCutMonitorHandler::fillCutMonitor(const MpdFemtoV0* v0, bool pass) {

  // Fill V0 cut monitors
  if (mCollectionsEmpty) return;

  MpdFemtoCutMonitorIterator iter;
  MpdFemtoBaseCutMonitor* CM;
  if (pass) {
    for (iter = mPassColl->begin(); iter != mPassColl->end(); iter++) {
      CM = *iter;
      CM->fill(v0);
    }
  } else {
    for (iter = mFailColl->begin(); iter != mFailColl->end(); iter++) {
      CM = *iter;
      CM->fill(v0);
    }
  }
}

//_________________
void MpdFemtoCutMonitorHandler::fillCutMonitor(const MpdFemtoKink* kink, bool pass) {

  if (mCollectionsEmpty) return;

  MpdFemtoCutMonitorIterator iter;
  MpdFemtoBaseCutMonitor* CM;

  if (pass) {
    for (iter = mPassColl->begin(); iter != mPassColl->end(); iter++) {
      CM = *iter;
      CM->fill(kink);
    }
  } else {
    for (iter = mFailColl->begin(); iter != mFailColl->end(); iter++) {
      CM = *iter;
      CM->fill(kink);
    }
  }
}

//_________________
void MpdFemtoCutMonitorHandler::fillCutMonitor(const MpdFemtoXi* xi, bool pass) {

  if (mCollectionsEmpty) return;

  MpdFemtoCutMonitorIterator iter;
  MpdFemtoBaseCutMonitor* CM;

  if (pass) {
    for (iter = mPassColl->begin(); iter != mPassColl->end(); iter++) {
      CM = *iter;
      CM->fill(xi);
    }
  } else {
    for (iter = mFailColl->begin(); iter != mFailColl->end(); iter++) {
      CM = *iter;
      CM->fill(xi);
    }
  }
}

//_________________
void MpdFemtoCutMonitorHandler::fillCutMonitor(const MpdFemtoPair* pair, bool pass) {

  if (mCollectionsEmpty) return;

  MpdFemtoCutMonitorIterator iter;
  MpdFemtoBaseCutMonitor* CM;

  if (pass) {
    for (iter = mPassColl->begin(); iter != mPassColl->end(); iter++) {
      CM = *iter;
      CM->fill(pair);
    }
  } else {
    for (iter = mFailColl->begin(); iter != mFailColl->end(); iter++) {
      CM = *iter;
      CM->fill(pair);
    }
  }
}

//_________________
void MpdFemtoCutMonitorHandler::fillCutMonitor(const MpdFemtoParticleCollection* partColl) {

  if (mCollectionsEmpty) return;

  MpdFemtoCutMonitorIterator iter;
  MpdFemtoBaseCutMonitor* CM;

  for (iter = mPassColl->begin(); iter != mPassColl->end(); iter++) {
    CM = *iter;
    CM->fill(partColl);
  }
}

//_________________
void MpdFemtoCutMonitorHandler::fillCutMonitor(const MpdFemtoEvent* event, const MpdFemtoParticleCollection* partColl) {

  // Fill event particle collection
  if (mCollectionsEmpty) return;

  MpdFemtoCutMonitorIterator iter;
  MpdFemtoBaseCutMonitor* CM;

  for (iter = mPassColl->begin(); iter != mPassColl->end(); iter++) {
    CM = *iter;
    CM->fill(event, partColl);
  }
}

//_________________
void MpdFemtoCutMonitorHandler::fillCutMonitor(const MpdFemtoParticleCollection *partColl1,
					       const MpdFemtoParticleCollection *partColl2) {
  // Fill event particle collection
  if (mCollectionsEmpty) return;

  MpdFemtoCutMonitorIterator iter;
  MpdFemtoBaseCutMonitor* CM;

  for (iter = mPassColl->begin(); iter != mPassColl->end(); iter++) {
    CM = *iter;
    CM->fill(partColl1, partColl2);
  }
}

//_________________
void MpdFemtoCutMonitorHandler::finish() {

  MpdFemtoCutMonitorIterator iter;
  for (iter = mPassColl->begin(); iter != mPassColl->end(); iter++) {
    (*iter)->finish();
  }
  for (iter = mFailColl->begin(); iter != mFailColl->end(); iter++) {
    (*iter)->finish();
  }
}

//_________________
void MpdFemtoCutMonitorHandler::addCutMonitor(MpdFemtoBaseCutMonitor* cutMoni1, MpdFemtoBaseCutMonitor* cutMoni2) {
  std::cout << "Adding cut monitor 1 ";
  std::cout << cutMoni1 << std::endl;
  std::cout << "Size of the first collection: " << mPassColl->size() << std::endl;
  mPassColl->push_back(cutMoni1);
  std::cout << "Adding cut monitor 2 ";
  std::cout << cutMoni1 << std::endl;
  mFailColl->push_back(cutMoni2);
  mCollectionsEmpty = false;
}

//_________________
void MpdFemtoCutMonitorHandler::addCutMonitor(MpdFemtoBaseCutMonitor* cutMoni) {
  //cout << " make a copy of the cutmonitor and push both into the collections " << endl;
  //cout << " not yet implemented" << endl;
  mPassColl->push_back(cutMoni);
  //cout << " only pass collection pushed" << endl;
  mCollectionsEmpty = false;
}

//_________________
void MpdFemtoCutMonitorHandler::addCutMonitorPass(MpdFemtoBaseCutMonitor* cutMoni) {
  mPassColl->push_back(cutMoni);
  mCollectionsEmpty = false;
}

//_________________
void MpdFemtoCutMonitorHandler::addCutMonitorFail(MpdFemtoBaseCutMonitor* cutMoni) {
  mFailColl->push_back(cutMoni);
  mCollectionsEmpty = false;
}

//_________________
MpdFemtoBaseCutMonitor* MpdFemtoCutMonitorHandler::passMonitor(int n) {

  if (static_cast<int> (mPassColl->size()) <= n) {
    return nullptr;
  }

  MpdFemtoCutMonitorIterator iter = mPassColl->begin();
  for (int i = 0; i < n; i++) {
    iter++;
  }
  return *iter;
}

//_________________
MpdFemtoBaseCutMonitor* MpdFemtoCutMonitorHandler::failMonitor(int n) {

  if (static_cast<int> (mFailColl->size()) <= n) {
    return nullptr;
  }

  MpdFemtoCutMonitorIterator iter = mFailColl->begin();
  for (int i = 0; i < n; i++) {
    iter++;
  }
  return *iter;
}

//_________________
TList *MpdFemtoCutMonitorHandler::getOutputList() {

  TList *tOutputList = new TList();

  for (auto &cut_monitor : *mPassColl) {
    TList *tLp = cut_monitor->getOutputList();

    TIter nextLp(tLp);
    while (TObject * obj = nextLp()) {
      tOutputList->Add(obj);
    }

    delete tLp;
  }

  for (auto &cut_monitor : *mFailColl) {
    TList *tLf = cut_monitor->getOutputList();

    TIter nextLf(tLf);
    while (TObject * obj = nextLf()) {
      tOutputList->Add(obj);
    }

    delete tLf;
  }

  return tOutputList;
}

//_________________
void MpdFemtoCutMonitorHandler::eventBegin(const MpdFemtoEvent* aEvent) {

  for (auto &cut_monitor : *mPassColl) {
    cut_monitor->eventBegin(aEvent);
  }

  for (auto &cut_monitor : *mFailColl) {
    cut_monitor->eventBegin(aEvent);
  }
}

//_________________
void MpdFemtoCutMonitorHandler::eventEnd(const MpdFemtoEvent* aEvent) {
  for (auto &cut_monitor : *mPassColl) {
    cut_monitor->eventEnd(aEvent);
  }

  for (auto &cut_monitor : *mFailColl) {
    cut_monitor->eventEnd(aEvent);
  }
}
