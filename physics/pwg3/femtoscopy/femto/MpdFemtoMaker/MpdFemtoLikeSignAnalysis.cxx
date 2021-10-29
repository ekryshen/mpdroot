//
// The class is a base class for like-sign analysis with z-binning
//

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseTrackCut.h"
#include "MpdFemtoBaseV0Cut.h"
// Infrastructure
#include "MpdFemtoLikeSignAnalysis.h"
#include "MpdFemtoParticleCollection.h"
#include "MpdFemtoPicoEventCollectionVector.h"
#include "MpdFemtoPicoEventCollectionVectorHideAway.h"

ClassImp(MpdFemtoLikeSignAnalysis)

/// This little function used to apply ParticleCuts (TrackCuts or V0Cuts)
/// and fill ParticleCollections of picoEvent it is called
/// from MpdFemtoAnalysis::ProcessEvent()
extern void fillHbtParticleCollection(MpdFemtoBaseParticleCut* partCut,
        MpdFemtoEvent* hbtEvent,
        MpdFemtoParticleCollection* partCollection, char doReshuffle);

//_________________
MpdFemtoLikeSignAnalysis::MpdFemtoLikeSignAnalysis(const unsigned int& bins, const double& min,
						   const double& max) : MpdFemtoAnalysis() {
  // Constructor
  mVertexBins = bins;
  mVertexZ[0] = min;
  mVertexZ[1] = max;
  mUnderFlow = 0;
  mOverFlow = 0;
  if (mMixingBuffer) delete mMixingBuffer;
  mPicoEventCollectionVectorHideAway =
    new MpdFemtoPicoEventCollectionVectorHideAway(mVertexBins, mVertexZ[0], mVertexZ[1]);
}

//_________________
MpdFemtoLikeSignAnalysis::MpdFemtoLikeSignAnalysis(const MpdFemtoLikeSignAnalysis& a) : MpdFemtoAnalysis(a) {
  // Copy constructor
  mVertexBins = a.mVertexBins;
  mVertexZ[0] = a.mVertexZ[0];
  mVertexZ[1] = a.mVertexZ[1];
  mUnderFlow = 0;
  mOverFlow = 0;
  if (mMixingBuffer) delete mMixingBuffer;
  mPicoEventCollectionVectorHideAway =
    new MpdFemtoPicoEventCollectionVectorHideAway(mVertexBins, mVertexZ[0], mVertexZ[1]);
}

//_________________
MpdFemtoLikeSignAnalysis& MpdFemtoLikeSignAnalysis::operator=(const MpdFemtoLikeSignAnalysis& a) {
  // Assignement operator
  if (this != &a) {
    mVertexBins = a.mVertexBins;
    mVertexZ[0] = a.mVertexZ[0];
    mVertexZ[1] = a.mVertexZ[1];
    mUnderFlow = 0;
    mOverFlow = 0;
    if (mMixingBuffer) delete mMixingBuffer;
    if (mPicoEventCollectionVectorHideAway) delete mPicoEventCollectionVectorHideAway;
    mPicoEventCollectionVectorHideAway =
      new MpdFemtoPicoEventCollectionVectorHideAway(mVertexBins, mVertexZ[0], mVertexZ[1]);
  }

  return *this;
}

//_________________
MpdFemtoLikeSignAnalysis::~MpdFemtoLikeSignAnalysis() {
  delete mPicoEventCollectionVectorHideAway;
  mPicoEventCollectionVectorHideAway = nullptr;
}

//_________________
MpdFemtoString MpdFemtoLikeSignAnalysis::report() {
  char Ctemp[200];
  std::cout << "MpdFemtoLikeSignAnalysis - constructing Report..." << std::endl;
  MpdFemtoString temp = "-----------\nHbt Analysis Report:\n";
  sprintf(Ctemp, "Events are mixed in %d bins in the range %E cm to %E cm.\n", mVertexBins, mVertexZ[0], mVertexZ[1]);
  temp += Ctemp;
  sprintf(Ctemp, "Events underflowing: %d\n", mUnderFlow);
  temp += Ctemp;
  sprintf(Ctemp, "Events overflowing: %d\n", mOverFlow);
  temp += Ctemp;
  sprintf(Ctemp, "Now adding MpdFemtoAnalysis(base) Report\n");
  temp += Ctemp;
  temp += "Adding MpdFemtoAnalysis(base) Report now:\n";
  temp += MpdFemtoAnalysis::report();
  temp += "-------------\n";
  MpdFemtoString returnThis = temp;
  return returnThis;
}

//_________________________
void MpdFemtoLikeSignAnalysis::processEvent(const MpdFemtoEvent* hbtEvent) {
  // Perform all the analysis tasks for a single event
  // get right mixing buffer
  double vertexZ = hbtEvent->primaryVertex().Z();
  mMixingBuffer = mPicoEventCollectionVectorHideAway->picoEventCollection(vertexZ);
  if (!mMixingBuffer) {
    if (vertexZ < mVertexZ[0]) mUnderFlow++;
    if (vertexZ > mVertexZ[1]) mOverFlow++;
    return;
  }

  // Startup for EbyE
  eventBegin(hbtEvent);

  // Event cut and event cut monitor
  bool tmpPassEvent = mEventCut->pass(hbtEvent);
  mEventCut->fillCutMonitor(hbtEvent, tmpPassEvent);
  if (tmpPassEvent) {
    mNeventsProcessed++;
    std::cout << "MpdFemtoLikeSignAnalysis::processEvent() - " << hbtEvent->trackCollection()->size();
    std::cout << " #track=" << hbtEvent->trackCollection()->size();
    // OK, analysis likes the event-- build a pico event from it,
    // using tracks the analysis likes...
    // This is what we will make pairs from and put in Mixing Buffer
    MpdFemtoPicoEvent* picoEvent = new MpdFemtoPicoEvent;

    fillHbtParticleCollection(mFirstParticleCut, (MpdFemtoEvent*) hbtEvent, picoEvent->firstParticleCollection(), mReshuffle);
    if (!(analyzeIdenticalParticles())) {
      fillHbtParticleCollection(mSecondParticleCut, (MpdFemtoEvent*) hbtEvent, picoEvent->secondParticleCollection(), mReshuffle);
    }

    std::cout << "   #particles in First, Second Collections: "
	      << picoEvent->firstParticleCollection()->size() << " "
	      << picoEvent->secondParticleCollection()->size() << std::endl;

    if (picoEvent->secondParticleCollection()->size() *
	picoEvent->firstParticleCollection()->size() == 0) {
      delete picoEvent;
      std::cout << "MpdFemtoLikeSignAnalysis - picoEvent deleted due to empty collection " << std::endl;
      return;
    }
    // OK, pico event is built
    // make real pairs...

    // Fabrice points out that we do not need to keep creating/deleting pairs all the time
    // We only ever need ONE pair, and we can just keep changing internal pointers
    // this should help speed things up
    MpdFemtoPair* ThePair = new MpdFemtoPair;

    MpdFemtoParticleIterator PartIter1;
    MpdFemtoParticleIterator PartIter2;
    MpdFemtoCorrFctnIterator CorrFctnIter;
    // Always
    MpdFemtoParticleIterator StartOuterLoop = picoEvent->firstParticleCollection()->begin();
    // Will be one less if identical
    MpdFemtoParticleIterator EndOuterLoop = picoEvent->firstParticleCollection()->end();
    MpdFemtoParticleIterator StartInnerLoop;
    MpdFemtoParticleIterator EndInnerLoop;
    // Only use First collection
    if (analyzeIdenticalParticles()) {
      // Iuter loop goes to next-to-last particle in First collection
      EndOuterLoop--;
      // Inner loop goes to last particle in First collection
      EndInnerLoop = picoEvent->firstParticleCollection()->end();
    } else {
      // Non-identical - loop over First and Second collections
      // Inner loop starts at first particle in Second collection
      StartInnerLoop = picoEvent->secondParticleCollection()->begin();
      // Inner loop goes to last particle in Second collection
      EndInnerLoop = picoEvent->secondParticleCollection()->end();
    }
    // real pairs
    for (PartIter1 = StartOuterLoop; PartIter1 != EndOuterLoop; PartIter1++) {
      if (analyzeIdenticalParticles()) {
	StartInnerLoop = PartIter1;
	StartInnerLoop++;
      }
      ThePair->setTrack1(*PartIter1);
      for (PartIter2 = StartInnerLoop; PartIter2 != EndInnerLoop; PartIter2++) {
	ThePair->setTrack2(*PartIter2);
	// The following lines have to be uncommented if you want pairCutMonitors
	// they are not in for speed reasons
	// bool tmpPassPair = mPairCut->Pass(ThePair);
	// mPairCut->FillCutMonitor(ThePair, tmpPassPair);
	// if ( tmpPassPair ) {
	if (mPairCut->pass(ThePair)) {
	  for (CorrFctnIter = mCorrFctnCollection->begin();
	       CorrFctnIter != mCorrFctnCollection->end(); CorrFctnIter++) {
	    MpdFemtoBaseLikeSignCorrFctn* CorrFctn = dynamic_cast<MpdFemtoBaseLikeSignCorrFctn*> (*CorrFctnIter);
	    if (CorrFctn) {
	      CorrFctn->addRealPair(ThePair);
	    }
	  }
	} //if ( mPairCut->pass(ThePair) )
      } //for (PartIter2 = StartInnerLoop; PartIter2!=EndInnerLoop; PartIter2++)
    } //for ( PartIter1=StartOuterLoop; PartIter1!=EndOuterLoop; PartIter1++ )

#ifdef STHBTDEBUG
    std::cout << "MpdFemtoLikeSignAnalysis::ProcessEvent() - reals done" << std::endl;
#endif

    MpdFemtoParticleIterator nextIter;
    MpdFemtoParticleIterator prevIter;

    // Like-sign first partilce collection pairs
    prevIter = EndOuterLoop;
    prevIter--;
    for (PartIter1 = StartOuterLoop; PartIter1 != prevIter; PartIter1++) {
      ThePair->setTrack1(*PartIter1);
      nextIter = PartIter1;
      nextIter++;
      for (PartIter2 = nextIter; PartIter2 != EndOuterLoop; PartIter2++) {
	ThePair->setTrack2(*PartIter2);
	// The following lines have to be uncommented if you want pairCutMonitors
	// they are not in for speed reasons
	// bool tmpPassPair = mPairCut->Pass(ThePair);
	// mPairCut->FillCutMonitor(ThePair, tmpPassPair);
	// if ( tmpPassPair ) {
	if (mPairCut->pass(ThePair)) {
	  for (CorrFctnIter = mCorrFctnCollection->begin();
	       CorrFctnIter != mCorrFctnCollection->end(); CorrFctnIter++) {
	    MpdFemtoBaseLikeSignCorrFctn* CorrFctn = dynamic_cast<MpdFemtoBaseLikeSignCorrFctn*> (*CorrFctnIter);
	    if (CorrFctn) {
	      CorrFctn->addLikeSignPositivePair(ThePair);
	    }
	  }
	} //if ( mPairCut->pass(ThePair) )
      } //for ( PartIter2 = nextIter; PartIter2!=EndOuterLoop; PartIter2++)
    } //for ( PartIter1=StartOuterLoop; PartIter1!=prevIter; PartIter1++)

#ifdef STHBTDEBUG
    std::cout << "MpdFemtoLikeSignAnalysis::processEvent() - like sign first collection done" << std::endl;
#endif

    // Like-sign second partilce collection pairs
    prevIter = EndInnerLoop;
    prevIter--;
    for (PartIter1 = StartInnerLoop; PartIter1 != prevIter; PartIter1++) {
      ThePair->setTrack1(*PartIter1);
      nextIter = PartIter1;
      nextIter++;
      for (PartIter2 = nextIter; PartIter2 != EndInnerLoop; PartIter2++) {
	ThePair->setTrack2(*PartIter2);
	// The following lines have to be uncommented if you want pairCutMonitors
	// they are not in for speed reasons
	// bool tmpPassPair = mPairCut->Pass(ThePair);
	// mPairCut->FillCutMonitor(ThePair, tmpPassPair);
	// if ( tmpPassPair ) {
	if (mPairCut->pass(ThePair)) {
	  for (CorrFctnIter = mCorrFctnCollection->begin();
	       CorrFctnIter != mCorrFctnCollection->end(); CorrFctnIter++) {
	    MpdFemtoBaseLikeSignCorrFctn* CorrFctn = dynamic_cast<MpdFemtoBaseLikeSignCorrFctn*> (*CorrFctnIter);
	    if (CorrFctn) {
	      CorrFctn->addLikeSignNegativePair(ThePair);
	    }
	  }
	} //if ( mPairCut->pass(ThePair) )
      } //for (PartIter2 = nextIter; PartIter2!=EndInnerLoop; PartIter2++)

      //for ( PartIter1=StartInnerLoop; PartIter1!=prevIter; PartIter1++)
#ifdef STHBTDEBUG
      std::cout << "MpdFemtoLikeSignAnalysis::processEvent() - like sign second collection done" << std::endl;
#endif

      if (mixingBufferFull()) {
#ifdef STHBTDEBUG
	std::cout << "Mixing Buffer is full - lets rock and roll" << std::endl;
#endif
      } else {
	std::cout << "Mixing Buffer not full -gotta wait " << mixingBuffer()->size() << std::endl;
      }

      if (mixingBufferFull()) {
	StartOuterLoop = picoEvent->firstParticleCollection()->begin();
	EndOuterLoop = picoEvent->firstParticleCollection()->end();
	MpdFemtoPicoEvent* storedEvent;
	MpdFemtoPicoEventIterator picoEventIter;
	for (picoEventIter = mixingBuffer()->begin();
	     picoEventIter != mixingBuffer()->end(); picoEventIter++) {
	  storedEvent = *picoEventIter;
	  if (analyzeIdenticalParticles()) {
	    StartInnerLoop = storedEvent->firstParticleCollection()->begin();
	    EndInnerLoop = storedEvent->firstParticleCollection()->end();
	  } else {
	    StartInnerLoop = storedEvent->secondParticleCollection()->begin();
	    EndInnerLoop = storedEvent->secondParticleCollection()->end();
	  }
	  for (PartIter1 = StartOuterLoop; PartIter1 != EndOuterLoop; PartIter1++) {
	    ThePair->setTrack1(*PartIter1);
	    for (PartIter2 = StartInnerLoop; PartIter2 != EndInnerLoop; PartIter2++) {
	      ThePair->setTrack2(*PartIter2);
	      // testing...	      cout << "ThePair defined... going to pair cut... ";
	      if (mPairCut->pass(ThePair)) {
		// testing...		cout << " ThePair passed PairCut... ";
		for (CorrFctnIter = mCorrFctnCollection->begin();
		     CorrFctnIter != mCorrFctnCollection->end(); CorrFctnIter++) {
		  MpdFemtoBaseLikeSignCorrFctn* CorrFctn = dynamic_cast<MpdFemtoBaseLikeSignCorrFctn*> (*CorrFctnIter);
		  if (CorrFctn) {
		    CorrFctn->addMixedPair(ThePair);
		    //cout << " ThePair has been added to MixedPair method " << endl;
		  } //if (CorrFctn)
		} // for ( CorrFctnIter=mCorrFctnCollection->begin();
	      } //if ( mPairCut->pass(ThePair) )
	    } //for (PartIter2=StartInnerLoop; PartIter2!=EndInnerLoop; PartIter2++)
	  } //for ( PartIter1=StartOuterLoop; PartIter1!=EndOuterLoop; PartIter1++)
	} //for ( picoEventIter=mixingBuffer()->begin(); picoEventIter!=mixingBuffer()->end(); picoEventIter++)

	// Now get rid of oldest stored pico-event in buffer.
	// This means (1) delete the event from memory, (2) "pop" the pointer to it from the MixingBuffer
	delete mixingBuffer()->back();
	mixingBuffer()->pop_back();
      } //if ( mixingBufferFull() )

      delete ThePair;
      // Store the current pico-event in buffer
      mixingBuffer()->push_front(picoEvent);
    } //if (tmpPassEvent)

    // Cleanup for EbyE
    eventEnd(hbtEvent);
    // cout << "MpdFemtoLikeSignAnalysis::ProcessEvent() - return to caller ... " << endl;
  }
}
