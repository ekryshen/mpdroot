//
// Femtoscopic analysis which mixes events with vertex z-binning
//

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseTrackCut.h"
#include "MpdFemtoBaseV0Cut.h"
#include "MpdFemtoBaseKinkCut.h"
// Infrastructure
#include "MpdFemtoVertexAnalysis.h"
#include "MpdFemtoParticleCollection.h"
#include "MpdFemtoPicoEventCollectionVector.h"
#include "MpdFemtoPicoEventCollectionVectorHideAway.h"

ClassImp(MpdFemtoVertexAnalysis)

//_________________
MpdFemtoVertexAnalysis::MpdFemtoVertexAnalysis(unsigned int bins, float zMin, float zMax) :
MpdFemtoAnalysis(),
mVertexBins(bins),
mOverFlow(0),
mUnderFlow(0) {
  // Parametrized constructor
  mVertexZ[0] = zMin;
  mVertexZ[1] = zMax;

  if (mMixingBuffer) {
    delete mMixingBuffer;
    mMixingBuffer = nullptr;
  }
  mPicoEventCollectionVectorHideAway =
    new MpdFemtoPicoEventCollectionVectorHideAway(mVertexBins, mVertexZ[0], mVertexZ[1]);
}

//_________________
MpdFemtoVertexAnalysis::MpdFemtoVertexAnalysis(const MpdFemtoVertexAnalysis& a) :
MpdFemtoAnalysis(a),
mVertexBins(a.mVertexBins),
mOverFlow(0),
mUnderFlow(0) {
  // Copy constructor
  mVertexZ[0] = a.mVertexZ[0];
  mVertexZ[1] = a.mVertexZ[1];

  if (mMixingBuffer) {
    delete mMixingBuffer;
    mMixingBuffer = nullptr;
  } // if (mMixingBuffer)

  mPicoEventCollectionVectorHideAway =
    new MpdFemtoPicoEventCollectionVectorHideAway(mVertexBins, mVertexZ[0], mVertexZ[1]);

  if (mVerbose) {
    std::cout << " MpdFemtoVertexAnalysis::MpdFemtoVertexAnalysis(const MpdFemtoVertexAnalysis& a) - analysis copied "
	      << std::endl;
  } // if (mVerbose)
}

//_________________
MpdFemtoVertexAnalysis& MpdFemtoVertexAnalysis::operator=(const MpdFemtoVertexAnalysis& a) {
  // Assignment operator
  if (this != &a) {
    MpdFemtoAnalysis::operator=(a);

    mVertexBins = a.mVertexBins;
    mVertexZ[0] = a.mVertexZ[0];
    mVertexZ[1] = a.mVertexZ[1];
    mUnderFlow = 0;
    mOverFlow = 0;

    if (mMixingBuffer) {
      delete mMixingBuffer;
      mMixingBuffer = nullptr;
    }
    mPicoEventCollectionVectorHideAway =
      new MpdFemtoPicoEventCollectionVectorHideAway(mVertexBins, mVertexZ[0], mVertexZ[1]);
  } // if ( this != &a )

  return *this;
}

//_________________
MpdFemtoVertexAnalysis::~MpdFemtoVertexAnalysis() {
  // Now delete every PicoEvent in the EventMixingBuffer and then the Buffer itself
  delete mPicoEventCollectionVectorHideAway;
}

//_________________
MpdFemtoString MpdFemtoVertexAnalysis::report() {
  // Prepare report fromt the execution
  if (mVerbose) {
    std::cout << "MpdFemtoVertexAnalysis - constructing report..." << std::endl;
  }

  TString report("-----------\nMpdFemtoVertexAnalysis report:\n");
  report += TString::Format("Events are mixed in %d bins in the range %E cm to %E cm.\n",
			    mVertexBins, mVertexZ[0], mVertexZ[1])
    + TString::Format("Events underflowing: %d\n", mUnderFlow)
    + TString::Format("Events overflowing: %d\n", mOverFlow)
    + TString::Format("Now adding MpdFemtoAnalysis(base) report\n");

  report += MpdFemtoAnalysis::report();

  return MpdFemtoString((const char *) report);
}

//_________________
void MpdFemtoVertexAnalysis::processEvent(const MpdFemtoEvent* hbtEvent) {

  if (mVerbose) {
    std::cout << " MpdFemtoVertexAnalysis::processEvent(const MpdFemtoEvent* hbtEvent) " << std::endl;
  }

  // Get right mixing buffer
  double vertexZ = hbtEvent->primaryVertex().Z();
  mMixingBuffer = mPicoEventCollectionVectorHideAway->picoEventCollection(vertexZ);

  if (!mMixingBuffer) {
    if (vertexZ < mVertexZ[0]) mUnderFlow++;
    if (vertexZ > mVertexZ[1]) mOverFlow++;
    return;
  }

  // Call processEvent() from MpdFemtoAnalysis-base
  MpdFemtoAnalysis::processEvent(hbtEvent);

  // NULL out the mixing buffer after event processed
  mMixingBuffer = nullptr;
}

//_________________
TList* MpdFemtoVertexAnalysis::listSettings() {
  TList *settings = MpdFemtoAnalysis::listSettings();
  settings->AddVector(new TObjString(TString::Format("MpdFemtoVertexAnalysis.vertex_z.bins=%d", mVertexBins)),
		      new TObjString(TString::Format("MpdFemtoVertexAnalysis.vertex_z.min=%f", mVertexZ[0])),
		      new TObjString(TString::Format("MpdFemtoVertexAnalysis.vertex_z.max=%f", mVertexZ[1])),
		      NULL);
  return settings;
}

//_________________________
TList* MpdFemtoVertexAnalysis::getOutputList() {
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
