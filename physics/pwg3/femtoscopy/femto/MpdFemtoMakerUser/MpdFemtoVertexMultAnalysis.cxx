//
// Perform femtoscopic analysis using z- and multiplicity binnings
//

// C++ headers
#include <iostream>

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseTrackCut.h"
#include "MpdFemtoBaseV0Cut.h"
#include "MpdFemtoBaseKinkCut.h"
#include "MpdFemtoBaseXiCut.h"
// Infrastructure
#include "MpdFemtoVertexMultAnalysis.h"
#include "MpdFemtoParticleCollection.h"
#include "MpdFemtoPicoEventCollectionVector.h"
#include "MpdFemtoPicoEventCollectionVectorHideAway.h"

// ROOT headers
#include "TString.h"

ClassImp(MpdFemtoVertexMultAnalysis)

//_________________
MpdFemtoVertexMultAnalysis::MpdFemtoVertexMultAnalysis(unsigned int binsVertex,
						       double minVertex,
						       double maxVertex,
						       unsigned int binsMult,
						       double minMult,
						       double maxMult) :
MpdFemtoAnalysis(),
mVertexZBins(binsVertex),
mOverFlowVertexZ(0),
mUnderFlowVertexZ(0),
mMultBins(binsMult),
mOverFlowMult(0),
mUnderFlowMult(0) {
  // Constructor
  mVertexZ[0] = minVertex;
  mVertexZ[1] = maxVertex;
  mMult[0] = minMult;
  mMult[1] = maxMult;

  // Print out warnings, will help user to fix these bugs
  if (minVertex >= maxVertex) {
    std::cout << "[WARNING] MpdFemtoVertexMultAnalysis - wrong z-vertex positions ("
	      << minVertex << " >= " << maxVertex << "). No events are expected to pass."
	      << std::endl;
  }

  if (minMult >= maxMult) {
    std::cout << "[WARNING] MpdFemtoVertexMultAnalysis - wrong multiplicity intervals ("
	      << minMult << " >= " << maxMult << "). No events are expected to pass."
	      << std::endl;
  }

  if (!mCorrFctnCollection) {
    mCorrFctnCollection = new MpdFemtoCorrFctnCollection;
  }

  // If the event collection was already create (it should NOT be) delete
  // before we allocate a new one
  if (mMixingBuffer) {
    delete mMixingBuffer;
    mMixingBuffer = nullptr;
  }

  mPicoEventCollectionVectorHideAway =
    new MpdFemtoPicoEventCollectionVectorHideAway(mVertexZBins, mVertexZ[0], mVertexZ[1],
						  mMultBins, mMult[0], mMult[1]);
}

//_________________
MpdFemtoVertexMultAnalysis::MpdFemtoVertexMultAnalysis(const MpdFemtoVertexMultAnalysis& a) :
MpdFemtoAnalysis(a),
mVertexZBins(a.mVertexZBins),
mOverFlowVertexZ(0),
mUnderFlowVertexZ(0),
mMultBins(a.mMultBins),
mOverFlowMult(0),
mUnderFlowMult(0) {
  // Copy constructor
  mVertexZ[0] = a.mVertexZ[0];
  mVertexZ[1] = a.mVertexZ[1];
  mMult[0] = a.mMult[0];
  mMult[1] = a.mMult[1];

  if (mMixingBuffer) {
    delete mMixingBuffer;
    mMixingBuffer = nullptr;
  }

  if (mPicoEventCollectionVectorHideAway) {
    delete mPicoEventCollectionVectorHideAway;
  }

  mPicoEventCollectionVectorHideAway =
    new MpdFemtoPicoEventCollectionVectorHideAway(mVertexZBins, mVertexZ[0], mVertexZ[1],
						  mMultBins, mMult[0], mMult[1]);

  if (mVerbose) {
    std::cout << "MpdFemtoVertexMultAnalysis::MpdFemtoVertexMultAnalysis(const MpdFemtoVertexMultAnalysis& a) - "
	      << "analysis copied" << std::endl;
  } // if (mVerbose)
}

//_________________
MpdFemtoVertexMultAnalysis& MpdFemtoVertexMultAnalysis::operator=(const MpdFemtoVertexMultAnalysis& a) {
  // Assignement operator
  if (this != &a) {

    // Allow parent class to copy the cuts and correlation functions
    MpdFemtoAnalysis::operator=(a);

    mVertexZBins = a.mVertexZBins;
    mMultBins = a.mMultBins;

    mVertexZ[0] = a.mVertexZ[0];
    mVertexZ[1] = a.mVertexZ[1];
    mUnderFlowVertexZ = 0;
    mOverFlowVertexZ = 0;

    mMult[0] = a.mMult[0];
    mMult[1] = a.mMult[1];
    mUnderFlowMult = 0;
    mOverFlowMult = 0;

    if (mMixingBuffer) {
      delete mMixingBuffer;
      mMixingBuffer = nullptr;
    }

    delete mPicoEventCollectionVectorHideAway;

    mPicoEventCollectionVectorHideAway =
      new MpdFemtoPicoEventCollectionVectorHideAway(mVertexZBins, mVertexZ[0], mVertexZ[1],
						    mMultBins, mMult[0], mMult[1]);
  } // if ( this != &a)

  return *this;
}

//_________________
MpdFemtoVertexMultAnalysis::~MpdFemtoVertexMultAnalysis() {
  // Destructor

  // Now delete every PicoEvent in the EventMixingBuffer and
  // then the Buffer itself
  delete mPicoEventCollectionVectorHideAway;
}

//_________________
MpdFemtoString MpdFemtoVertexMultAnalysis::report() {
  // Prepare a report of the execution
  if (mVerbose) {
    std::cout << "MpdFemtoVertexMultAnalysis - constructing report..." << std::endl;
  }

  TString report("-----------\nHbt   MpdFemtoVertexMultAnalysis Report:\n");

  report += TString::Format("Events are mixed in %d VertexZ bins in the range %E cm to %E cm.\n",
			    mVertexZBins, mVertexZ[0], mVertexZ[1])
    + TString::Format("Events underflowing: %d\n", mUnderFlowVertexZ)
    + TString::Format("Events overflowing: %d\n", mOverFlowVertexZ)
    + TString::Format("Events are mixed in %d Mult bins in the range %E cm to %E cm.\n",
		      mMultBins, mMult[0], mMult[1])
    + TString::Format("Events underflowing: %d\n", mUnderFlowMult)
    + TString::Format("Events overflowing: %d\n", mOverFlowMult)
    + TString::Format("Now adding MpdFemtoAnalysis(base) report\n")
    + MpdFemtoAnalysis::report();

  return MpdFemtoString((const char *) report);
}

//_________________
TList* MpdFemtoVertexMultAnalysis::listSettings() {
  TList *settings = MpdFemtoAnalysis::listSettings();
  settings->AddVector(new TObjString(TString::Format("MpdFemtoVertexMultAnalysis.vertex_z.bins=%d", mVertexZBins)),
		      new TObjString(TString::Format("MpdFemtoVertexMultAnalysis.vertex_z.min=%f", mVertexZ[0])),
		      new TObjString(TString::Format("MpdFemtoVertexMultAnalysis.vertex_z.max=%f", mVertexZ[1])),
		      new TObjString(TString::Format("MpdFemtoVertexMultAnalysis.multiplicity.bins=%d", mMultBins)),
		      new TObjString(TString::Format("MpdFemtoVertexMultAnalysis.multiplicity.min=%f", mMult[0])),
		      new TObjString(TString::Format("MpdFemtoVertexMultAnalysis.multiplicity.max=%f", mMult[1])),
		      NULL);
  return settings;
}

//_________________________
TList* MpdFemtoVertexMultAnalysis::getOutputList() {
  // Collect the list of output objects to be written
  TList *tOutputList = new TList();

  TList *p1Cut = mFirstParticleCut->getOutputList();
  TListIter nextp1(p1Cut);
  while (TObject * obj = nextp1.Next()) {
    tOutputList->Add(obj);
  }
  delete p1Cut;

  if (mSecondParticleCut != mFirstParticleCut) {
    TList *p2Cut = mSecondParticleCut->getOutputList();

    TIter nextp2(p2Cut);
    while (TObject * obj = nextp2()) {
      tOutputList->Add(obj);
    }
    delete p2Cut;
  } //if (fSecondParticleCut != fFirstParticleCut)

  TList *pairCut = mPairCut->getOutputList();
  TIter nextpair(pairCut);
  while (TObject * obj = nextpair()) {
    tOutputList->Add(obj);
  }
  delete pairCut;

  TList *eventCut = mEventCut->getOutputList();
  TIter nextevent(eventCut);
  while (TObject * obj = nextevent()) {
    tOutputList->Add(obj);
  }
  delete eventCut;

  for (auto &cf : *mCorrFctnCollection) {
    TList *tListCf = cf->getOutputList();

    TIter nextListCf(tListCf);
    while (TObject * obj = nextListCf()) {
      tOutputList->Add(obj);
    }
    delete tListCf;
  }

  return tOutputList;
}

//_________________
void MpdFemtoVertexMultAnalysis::processEvent(const MpdFemtoEvent* hbtEvent) {
  // Perform event processing in bins of z vertex and multiplicity

  // Find the correct mixing buffer
  double vertexZ = hbtEvent->primaryVertex().Z();
  double mult = hbtEvent->refMult();
  mMixingBuffer = mPicoEventCollectionVectorHideAway->picoEventCollection(vertexZ, mult);
  if (!mMixingBuffer) {
    if (vertexZ < mVertexZ[0]) mUnderFlowVertexZ++;
    if (vertexZ > mVertexZ[1]) mOverFlowVertexZ++;
    if (mult < mMult[0]) mUnderFlowMult++;
    if (mult > mMult[1]) mOverFlowMult++;
    return;
  }

  // Call ProcessEvent() from MpdFemtoAnalysis-base
  MpdFemtoAnalysis::processEvent(hbtEvent);

  // NULL out the mixing buffer after event processed
  mMixingBuffer = nullptr;
}
