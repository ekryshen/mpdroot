//
// An example of the most basic (concrete) analysis
//

// C++ headers
#include <string>
#include <iostream>
#include <iterator>
#include <numeric>
#include <algorithm>
#include <random>

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseTrackCut.h"
#include "MpdFemtoBaseV0Cut.h"
#include "MpdFemtoBaseXiCut.h"
#include "MpdFemtoBaseKinkCut.h"
// Infrastructure
#include "MpdFemtoAnalysis.h"

// ROOT headers
#include "TObject.h"
#include "TString.h"

ClassImp(MpdFemtoAnalysis)

MpdFemtoBaseEventCut* copyTheCut(MpdFemtoBaseEventCut*);
MpdFemtoBaseParticleCut* copyTheCut(MpdFemtoBaseParticleCut*);
MpdFemtoBasePairCut* copyTheCut(MpdFemtoBasePairCut*);
MpdFemtoBaseCorrFctn* copyTheCorrFctn(MpdFemtoBaseCorrFctn*);

/// Generalized particle collection filler function - called by
/// FillParticleCollection()
///
/// This function loops over the track_collection, calling the cut's Pass
/// method on each track. If it passes, a new MpdFemtoParticle is constructed
/// from the track's pointer (or more generally, the item returned by
/// dereferencing the container's iterator) and the cut's expected mass.
///
/// This templated function accepts a track cut, track collection, and an
/// MpdFemtoParticleCollection (which points to the output) as input. The types
/// of the tracks are determined by the template paramters, which should be
/// automatically detected by argument inspection (you don't need to specify).
///
/// The original function also specified the iterator type, but since all
/// containers are standard STL containers, it now infers the type as
/// TrackCollectionType::iterator. If the track collections change to some
/// other type, it is recommended to add TrackCollectionIterType to the
/// template list, and add the appropriate type to the function calls in
/// FillParticleCollection.
template <class TrackCollectionType, class TrackCutType>
void doFillParticleCollection(TrackCutType *cut,
			      TrackCollectionType *track_collection,
			      MpdFemtoParticleCollection *output, char doReshuffle) {
  MpdFemtoParticleCollection* initPartCollection = new MpdFemtoParticleCollection();
  for (const auto &track : *track_collection) {
    const Bool_t track_passes = cut->pass(track);
    cut->fillCutMonitor(track, track_passes);
    if (track_passes) {
      if (doReshuffle == 2)
	initPartCollection->push_back(new MpdFemtoParticle(track, cut->mass()));
      else
	output->push_back(new MpdFemtoParticle(track, cut->mass()));

    } //if (track_passes)
  } //for (const auto &track : *track_collection)

    //reshuffling implementation by means of reshuffling initPartCollection iterators in vector iterColl
  if (initPartCollection && doReshuffle == 2) {
    unsigned int collectionSize = initPartCollection->size();
    std::vector<MpdFemtoParticleIterator>* iterColl = new std::vector<MpdFemtoParticleIterator>(collectionSize);
    std::iota(iterColl->begin(), iterColl->end(), initPartCollection->begin()); //fill iterColl by iterators

    //  debug---->
    //    unsigned int iterCollSize = iterColl->size();
    //    if(collectionSize!=iterCollSize){
    //      std::cout<<"[ERROR]: Size of particles collection and size of iterators for this collection are different"<<
    //				"\nSomthing is wrong with function std::iota ..."<<std::endl;
    //      initPartCollection->clear();
    //      return;
    //    }
    //  <----debug

    // Random generator initialization
    std::random_device rd;
    // Random generator initialization
    std::mt19937 gen(rd());
    std::shuffle(iterColl->begin(), iterColl->end(), gen);

    //  debug---->
    //    iterCollSize = iterColl->size();
    //    if(collectionSize!=iterCollSize){
    //      std::cout<<"[ERROR]: Size of iterators vector is changed after reshuffling"<<
    //				"\nSomthing is wrong with function std::reshuffle..."<<std::endl;
    //      initPartCollection->clear();
    //      return;
    //    }
    //  <----debug

    for (vector <MpdFemtoParticleIterator>::iterator i = iterColl->begin(); i != iterColl->end(); i++) {
      output->push_back((*(*i)));
    }

    //  debug---->
    //  unsigned int outputSize = output->size();
    //  if(collectionSize!=outputSize){
    //    std::cout<<"[ERROR]: Size of particles collection is changed\n"<<
    //               "during convertion from iterators vector to result particles collection"<<std::endl;
    //    initPartCollection->clear();
    //    return;
    //  }
    //  <----debug

    // debug printing BEGIN
    // Let's check we don't lose any particle because of reshuffling process
    // WARNING: Plain and really slow realization! Uncomment for debugging only
    // Idea: at first we verify that resulting collection contains all poiners MpdFemtoParticle*
    // from initial collection.
    // If it's right reshuffling process is correct, if it's not we should check vector of
    // iterators that participate 
    // in reshuffling: 1. vector should be filled correctly; 2. all iterators are in this
    // vector after reshuffling algorithm.

    //  //Sizes of all collections have been already checked in code wrapped by "//debug--->	//<---debug"
    //  bool found=0;
    //  unsigned int k=0, l=0;
    //  std::cout<<"Reshuffling result:"<<std::endl;
    //  for(MpdFemtoParticleIterator i = initPartCollection->begin(); i != initPartCollection->end(); i++){
    //    found=0;
    //    l=0;
    //    for(MpdFemtoParticleIterator j = output->begin(); j != output->end(); j++){
    //      if((*i)==(*j)){
    //        std::cout<<k<<"-"<<l<<"		";
    //        found=1;
    //        break;
    //      }
    //      l++;
    //    }//for(MpdFemtoParticleIterator j = output->begin(); j != output->end(); j++)
    //    if(found==0){
    //      std::cout<<"[ERROR]: "<<k<<" particle can't be found. Let's check inner processes...'"<<std::endl;
    //      break;
    //    }
    //    k++;
    //  }//for(MpdFemtoParticleIterator i = initPartCollection->begin(); i != initPartCollection->end(); i++)
    //  std::cout<<std::endl;

    //  if(found==0){
    //    //Vector of iterators verification
    //    k=0;
    //    for(MpdFemtoParticleIterator i = initPartCollection->begin(); i != initPartCollection->end(); i++){
    //      if((*i) != *(iterColl->at(k)) ){
    //        std::cout<<"[ERROR]: Wrong filling of iterators vector. Broken at "<<k<<" particle"<<
    //                   "\nSomthing is wrong with function std::iota ..."<<std::endl;
    //        initPartCollection->clear();
    //        return;          
    //      }
    //      k++;      
    //    }//for(unsigned int i=0; i<collectionSize;i++)

    //   //If this point is reached it means some problems are in std::shuffle.
    //   std::cout<<"[ERROR]: std::shuffle loses some particles.\nYou lose."<<std::endl;
    //  }//if(found==0) 
    //  else{
    //    std::cout<<"You win!\nReshuffling process is ok.\n"<<std::endl;
    //  }  
    // //debug printing END

    initPartCollection->clear();
  }//if(initPartCollection && doReshuffle==2)
}

/// This little function is used to apply ParticleCuts (TrackCuts or V0Cuts) and
/// fill ParticleCollections from tacks in picoEvent. It is called from
/// MpdFemtoAnalysis::processEvent().
///
/// The actual loop implementation has been moved to the collection-generic
/// doFillParticleCollection() function
void fillHbtParticleCollection(MpdFemtoBaseParticleCut* partCut,
			       MpdFemtoEvent* hbtEvent,
			       MpdFemtoParticleCollection* partCollection, char doReshuffle) {
  // Selection of the particle types: Track, V0, Kink
  switch (partCut->type()) {

  case hbtTrack:
    // Cut is cutting on Tracks
    doFillParticleCollection((MpdFemtoBaseTrackCut*) partCut,
			     hbtEvent->trackCollection(),
			     partCollection, doReshuffle);
    break;

  case hbtV0:
    // Cut is cutting on V0s
    doFillParticleCollection((MpdFemtoBaseV0Cut*) partCut,
			     hbtEvent->v0Collection(),
			     partCollection, doReshuffle);
    break;

  case hbtXi:
    // Cut is cutting on Xis
    doFillParticleCollection((MpdFemtoBaseXiCut*) partCut,
			     hbtEvent->xiCollection(),
			     partCollection, doReshuffle);
    break;

  case hbtKink:
    // Cut is cutting on Kinks
    doFillParticleCollection((MpdFemtoBaseKinkCut*) partCut,
			     hbtEvent->kinkCollection(),
			     partCollection, doReshuffle);
    break;
  default:
    std::cout << "fillHbtParticleCollection function (in MpdFemtoAnalysis.cxx): "
	      << "Undefined Particle Cut type!!!" << partCut->type() << std::endl;
  } //switch (partCut->Type())

  partCut->fillCutMonitor(hbtEvent, partCollection);
}

//_________________
MpdFemtoAnalysis::MpdFemtoAnalysis() :
  mPicoEventCollectionVectorHideAway(nullptr),
  mPairCut(nullptr), mCorrFctnCollection(nullptr),
  mEventCut(nullptr), mFirstParticleCut(nullptr),
  mSecondParticleCut(nullptr), mMixingBuffer(nullptr),
  mPicoEvent(nullptr), mNumEventsToMix(0), mNeventsProcessed(0),
  mMinSizePartCollection(0), mReshuffle(1),
  mVerbose(false) {
  // Default constructor
  mCorrFctnCollection = new MpdFemtoCorrFctnCollection;
  mMixingBuffer = new MpdFemtoPicoEventCollection;
}

//_________________
MpdFemtoAnalysis::MpdFemtoAnalysis(const MpdFemtoAnalysis& a) :
  MpdFemtoBaseAnalysis(),
  mPicoEventCollectionVectorHideAway(nullptr),
  mPairCut(nullptr),
  mCorrFctnCollection(nullptr),
  mEventCut(nullptr),
  mFirstParticleCut(nullptr),
  mSecondParticleCut(nullptr),
  mMixingBuffer(nullptr),
  mPicoEvent(nullptr),
  mNumEventsToMix(a.mNumEventsToMix),
  mNeventsProcessed(0),
  mMinSizePartCollection(a.mMinSizePartCollection),
  mReshuffle(a.mReshuffle),
  mVerbose(a.mVerbose) {
  // Copy constructor
  const char msg_template[] = " MpdFemtoAnalysis::MpdFemtoAnalysis(const MpdFemtoAnalysis& a) - %s";
  const char warn_template[] = " [WARNING] MpdFemtoAnalysis::MpdFemtoAnalysis(const MpdFemtoAnalysis& a)] %s";

  mCorrFctnCollection = new MpdFemtoCorrFctnCollection;
  mMixingBuffer = new MpdFemtoPicoEventCollection;

  // Clone the event cut
  mEventCut = a.mEventCut->clone();
  if (mEventCut) {
    setEventCut(mEventCut);
    if (mVerbose) {
      std::cout << TString::Format(msg_template, "Event cut set") << std::endl;
    } // if( mVerbose )
  } else {
    std::cerr << Form(warn_template, "Cannot clone event cut!") << std::endl;
  }

  // Clone the first particle cut
  mFirstParticleCut = a.mFirstParticleCut->clone();
  if (mFirstParticleCut) {
    setFirstParticleCut(mFirstParticleCut);
    if (mVerbose) {
      std::cout << TString::Format(msg_template, "First particle cut set") << std::endl;
    } // if( mVerbose )
  } else {
    std::cerr << TString::Format(warn_template, "Cannot clone first particle cut!") << std::endl;
  }

  // Clone the second particle cut (if it exists)
  mSecondParticleCut = ((a.mFirstParticleCut == a.mSecondParticleCut) ?
			mFirstParticleCut : a.mSecondParticleCut);
  if (mSecondParticleCut) {
    setSecondParticleCut(mSecondParticleCut);
    if (mVerbose) {
      std::cout << TString::Format(msg_template, "Second particle cut set") << std::endl;
    } // if( mVerbose )
  } else {
    std::cerr << TString::Format(warn_template, "Cannot clone second particle cut!") << std::endl;
  }

  // Clone the pair cut
  mPairCut = a.mPairCut->clone();
  if (mPairCut) {
    setPairCut(mPairCut);
    if (mVerbose) {
      std::cout << TString::Format(msg_template, "Pair cut set") << std::endl;
    } // if( mVerbose )
  } else {
    std::cerr << TString::Format(warn_template, "Cannot clone pair cut!") << std::endl;
  }

  MpdFemtoCorrFctnIterator iter;
  for (iter = a.mCorrFctnCollection->begin();
       iter != a.mCorrFctnCollection->end(); iter++) {

    if (mVerbose) {
      std::cout << TString::Format(msg_template, "looking for correlation functions") << std::endl;
    } // if( mVerbose )
    MpdFemtoBaseCorrFctn* fctn = (*iter)->clone();
    if (fctn) {
      addCorrFctn(fctn);
    } else {
      std::cout << TString::Format(msg_template, "correlation function not found") << std::endl;
    }
  }

  if (mVerbose) {
    std::cout << TString::Format(msg_template, "Analysis copied") << std::endl;
  }
}

//_________________
MpdFemtoAnalysis& MpdFemtoAnalysis::operator=(const MpdFemtoAnalysis& ana) {

  // Assignment operator
  if (this != &ana) {

    const char warn_template[] = "MpdFemtoAnalysis& MpdFemtoAnalysis::operator=(const MpdFemtoAnalysis& ana) : %s";

    // Clear second particle cut to avoid double delete
    if (mFirstParticleCut == mSecondParticleCut) {
      mSecondParticleCut = nullptr;
    }
    if (mEventCut) delete mEventCut;
    if (mPairCut) delete mPairCut;
    if (mFirstParticleCut) delete mFirstParticleCut;
    if (mSecondParticleCut) delete mSecondParticleCut;

    // Clear correlation functions
    if (mCorrFctnCollection) {
      for (MpdFemtoCorrFctnIterator iIter = mCorrFctnCollection->begin();
	   iIter != mCorrFctnCollection->end(); iIter++) {
	delete *iIter;
      }
      mCorrFctnCollection->clear();
    }//if (mCorrFctnCollection)
    else {
      std::cerr << TString::Format(warn_template, "mCorrFctnCollection is a nullptr") << std::endl;
      mCorrFctnCollection = new MpdFemtoCorrFctnCollection;
    }

    // Clear mixing buffer
    if (mMixingBuffer) {
      for (MpdFemtoPicoEventIterator iIter = mMixingBuffer->begin();
	   iIter != mMixingBuffer->end(); iIter++) {
	delete *iIter;
      }
      mMixingBuffer->clear();
    }//if (mMixingBuffer)
    else {
      std::cerr << TString::Format(warn_template, "mMixingBuffer is a nullptr") << std::endl;
      mMixingBuffer = new MpdFemtoPicoEventCollection;
    }

    // Copy objects
    mPairCut = ana.mPairCut->clone();
    mEventCut = ana.mEventCut->clone();
    mFirstParticleCut = ana.mFirstParticleCut->clone();
    mSecondParticleCut = ((ana.mFirstParticleCut == ana.mSecondParticleCut) ?
			  mFirstParticleCut : ana.mSecondParticleCut->clone());

    // Check that cuts have just been coppied
    if (mEventCut) {
      setEventCut(mEventCut);
    } else {
      std::cerr << TString::Format(warn_template, "Could not copy the event cut") << std::endl;
    }

    if (mPairCut) {
      setPairCut(mPairCut);
    } else {
      std::cerr << TString::Format(warn_template, "Could not copy the pair cut") << std::endl;
    }

    if (mFirstParticleCut) {
      setFirstParticleCut(mFirstParticleCut);
    } else {
      std::cerr << TString::Format(warn_template, "Could not copy the first particle cut") << std::endl;
    }

    if (mSecondParticleCut) {
      setSecondParticleCut(mSecondParticleCut);
    } else {
      std::cerr << TString::Format(warn_template, "Could not copy the second particle cut") << std::endl;
    }

    // Copy correlation functions
    for (auto &cf : *ana.mCorrFctnCollection) {
      MpdFemtoBaseCorrFctn *fctn = cf->clone();
      if (fctn) {
	addCorrFctn(fctn);
      }
    }

    // Copy some numbers
    mNumEventsToMix = ana.mNumEventsToMix;
    mMinSizePartCollection = ana.mMinSizePartCollection;
    mReshuffle = ana.mReshuffle;
    mVerbose = ana.mVerbose;
  } //if ( this != &ana )

  return *this;
}

//_________________
MpdFemtoAnalysis::~MpdFemtoAnalysis() {

  std::cout << "MpdFemtoAnalysis::~MpdFemtoAnalysis()" << std::endl;

  // Delete event cut
  if (mEventCut) delete mEventCut;
  mEventCut = nullptr;
  // Double-delete protection
  if (mFirstParticleCut == mSecondParticleCut) {
    mSecondParticleCut = nullptr;
  }
  if (mFirstParticleCut) delete mFirstParticleCut;
  mFirstParticleCut = nullptr;
  if (mSecondParticleCut) delete mSecondParticleCut;
  mSecondParticleCut = nullptr;
  if (mPairCut) delete mPairCut;
  mPairCut = nullptr;

  // Now delete every CorrFunction in the Collection, and then the Collection itself
  if (mCorrFctnCollection) {
    MpdFemtoCorrFctnIterator iter;
    for (iter = mCorrFctnCollection->begin();
	 iter != mCorrFctnCollection->end(); iter++) {
      delete *iter;
    }
    delete mCorrFctnCollection;
  } //if( mCorrFctnCollection )

    // Now delete every PicoEvent in the EventMixingBuffer and then the Buffer itself
  if (mMixingBuffer) {
    MpdFemtoPicoEventIterator piter;
    for (piter = mMixingBuffer->begin();
	 piter != mMixingBuffer->end(); piter++) {
      delete *piter;
    }
    delete mMixingBuffer;
  } //if (mMixingBuffer)
}

//_________________
MpdFemtoBaseCorrFctn* MpdFemtoAnalysis::corrFctn(int n) {

  // Return pointer to n-th correlation function

  if ((n < 0) || (static_cast<size_t> (n) > mCorrFctnCollection->size())) {
    return nullptr;
  }
  MpdFemtoCorrFctnIterator iter = mCorrFctnCollection->begin();
  for (int i = 0; i < n; i++) {
    iter++;
  }
  return *iter;
}

//_________________________
MpdFemtoString MpdFemtoAnalysis::reshReport() {
  string temp = "\nParticles reshuffling: ";
  switch (MpdFemtoAnalysis::reshuffle()) {
  case 0:
    {
      temp += "OFF";
      break;
    }
  case 1:
    {
      temp += "ON, swapping particles in pairs\n";
      break;
    }
  case 2:
    {
      temp += "ON, Fisher-Yates algorithm for particle collection\n";
      break;
    }
  default:
    temp += "[ERROR] MpdFemtoAnalsysis::mReshuffle has invalid value\n";
  }
  MpdFemtoString returnThis = temp;
  return returnThis;
}

//_________________
MpdFemtoString MpdFemtoAnalysis::report() {

  // Create a simple report from the analysis execution
  std::cout << "MpdFemtoAnalysis - constructing Report..."
            << std::endl;
  string temp = "-----------\nHbt Analysis Report:\n";
  temp += reshReport();
  temp += "\nEvent Cuts:\n";
  temp += mEventCut->report();
  temp += "\nParticle Cuts - First Particle:\n";
  temp += mFirstParticleCut->report();
  temp += "\nParticle Cuts - Second Particle:\n";
  temp += mSecondParticleCut->report();
  temp += "\nPair Cuts:\n";
  temp += mPairCut->report();
  temp += "\nCorrelation Functions:\n";
  MpdFemtoCorrFctnIterator iter;
  if (mCorrFctnCollection->size() == 0) {
    std::cout << "MpdFemtoAnalysis-Warning : no correlations functions in this analysis "
	      << std::endl;
  }
  for (iter = mCorrFctnCollection->begin(); iter != mCorrFctnCollection->end(); iter++) {
    temp += (*iter)->report();
    temp += "\n";
  }
  temp += "-------------\n";
  MpdFemtoString returnThis = temp;
  return returnThis;
}

//_________________
void MpdFemtoAnalysis::processEvent(const MpdFemtoEvent* hbtEvent) {

  // Add event to processed events

  // We will get a new pico event, if not prevent
  // correlation function to access old pico event
  mPicoEvent = nullptr;
  addEventProcessed();

  // Startup for EbyE
  eventBegin(hbtEvent);

  // Event cut and event cut monitor
  bool tmpPassEvent = mEventCut->pass(hbtEvent);
  if (!tmpPassEvent) {
    mEventCut->fillCutMonitor(hbtEvent, tmpPassEvent);
    // Cleanup for EbyE
    eventEnd(hbtEvent);
    return;
  }
  // Analysis likes the event -- build a pico event from it, using tracks the
  // analysis likes. This is what we will make pairs from and put in Mixing
  // Buffer.
  // No memory leak. We will delete picoevent when they come out
  // of the mixing buffer
  mPicoEvent = new MpdFemtoPicoEvent;

  MpdFemtoParticleCollection *collection1 = mPicoEvent->firstParticleCollection();
  MpdFemtoParticleCollection *collection2 = mPicoEvent->secondParticleCollection();
  // Sanity check that both collections exist
  if (collection1 == nullptr || collection2 == nullptr) {
    std::cout << "MpdFemtoAnalysis::processEvent - new PicoEvent is missing particle collections!"
	      << std::endl;
    eventEnd(hbtEvent);
    delete mPicoEvent;
    return;
  }
  // Subroutine fills fPicoEvent'a FirstParticleCollection with tracks from
  // hbtEvent which pass mFirstParticleCut. Uses cut's "type()" to determine
  // which track collection to pull from hbtEvent.

  fillHbtParticleCollection(mFirstParticleCut, (MpdFemtoEvent*) hbtEvent,
			    mPicoEvent->firstParticleCollection(), mReshuffle);
  // In case of non-identical particles
  if (!(analyzeIdenticalParticles())) {
    fillHbtParticleCollection(mSecondParticleCut, (MpdFemtoEvent*) hbtEvent,
			      mPicoEvent->secondParticleCollection(), mReshuffle);
  }

  if (mVerbose) {
    std::cout << " MpdFemtoAnalysis::processEvent - #particles in First, Second Collections: "
	      << mPicoEvent->firstParticleCollection()->size() << " "
	      << mPicoEvent->secondParticleCollection()->size() << std::endl;
  } // if ( mVerbose )

  const unsigned int coll_1_size = collection1->size();
  const unsigned int coll_2_size = collection2->size();
  const bool coll_1_size_passes = coll_1_size >= mMinSizePartCollection;
  const bool coll_2_size_passes = analyzeIdenticalParticles() || (coll_2_size >= mMinSizePartCollection);

  tmpPassEvent = (tmpPassEvent && coll_1_size_passes && coll_2_size_passes);

  // Fill event cut monitor
  mEventCut->fillCutMonitor(hbtEvent, tmpPassEvent);

  // Stop here if event did not pass cuts
  if (!tmpPassEvent) {
    eventEnd(hbtEvent);
    delete mPicoEvent;
    return;
  }

  //------ Make real pairs. If identical, make pairs for one collection ------//
  if (analyzeIdenticalParticles()) {
    makePairs("real", mPicoEvent->firstParticleCollection());
  } else {
    makePairs("real", mPicoEvent->firstParticleCollection(),
	      mPicoEvent->secondParticleCollection());
  }

  if (mVerbose) {
    std::cout << "MpdFemtoAnalysis::processEvent() - reals done ";
  }

  //---- Make pairs for mixed events, looping over events in mixingBuffer ----//

  MpdFemtoPicoEvent* storedEvent;
  MpdFemtoPicoEventIterator mPicoEventIter;
  for (mPicoEventIter = mixingBuffer()->begin();
       mPicoEventIter != mixingBuffer()->end(); mPicoEventIter++) {

    storedEvent = *mPicoEventIter;

    if (analyzeIdenticalParticles()) {
      makePairs("mixed", mPicoEvent->firstParticleCollection(),
		storedEvent->firstParticleCollection());
    } else {
      makePairs("mixed", mPicoEvent->firstParticleCollection(),
		storedEvent->secondParticleCollection());

      makePairs("mixed", storedEvent->firstParticleCollection(),
		mPicoEvent->secondParticleCollection());
    }
  } //for ( mPicoEventIter=MixingBuffer()->begin();	mPicoEventIter!=MixingBuffer()->end(); mPicoEventIter++)

  if (mVerbose) {
    std::cout << " - mixed done   " << std::endl;
  } // if (mVerbose)

    //--------- If mixing buffer is full, delete oldest event ---------//
  if (mixingBufferFull()) {
    delete mixingBuffer()->back();
    mixingBuffer()->pop_back();
  }

  //-------- Add current event (mPicoEvent) to mixing buffer --------//
  mixingBuffer()->push_front(mPicoEvent);

  // Cleanup for EbyE
  eventEnd(hbtEvent);
  //std::cout << "MpdFemtoAnalysis::ProcessEvent() - return to caller ... " << std::endl;
}

//_________________________
void MpdFemtoAnalysis::makePairs(const char* typeIn,
				 MpdFemtoParticleCollection *partCollection1,
				 MpdFemtoParticleCollection *partCollection2) {

  // Build pairs, check pair cuts, and call CFs' AddRealPair() or
  // AddMixedPair() methods. If no second particle collection is
  // specfied, make pairs within first particle collection.
  const string type = typeIn;

  bool swpart = mNeventsProcessed % 2;

  // Create the pair outside  the loop
  MpdFemtoPair* ThePair = new MpdFemtoPair;

  MpdFemtoCorrFctnIterator CorrFctnIter;
  MpdFemtoParticleIterator PartIter1, PartIter2;

  // Setup iterator ranges
  //
  // The outer loop alway starts at beginning of particle collection 1.
  // * If we are iterating over both particle collections, then the loop simply
  // runs through both from beginning to end.
  // * If we are only iterating over one particle collection, the inner loop
  // loops over all particles between the outer iterator and the end of the
  // collection. The outer loop must skip the last entry of the list.
  MpdFemtoParticleIterator StartOuterLoop = partCollection1->begin();
  MpdFemtoParticleIterator EndOuterLoop = partCollection1->end();
  MpdFemtoParticleIterator StartInnerLoop;
  MpdFemtoParticleIterator EndInnerLoop;

  if (partCollection2) { // Two collections:
    StartInnerLoop = partCollection2->begin(); //   Full inner & outer loops
    EndInnerLoop = partCollection2->end(); //
  } else { // One collection:
    EndOuterLoop--; //   Outer loop goes to next-to-last particle
    EndInnerLoop = partCollection1->end(); //   Inner loop goes to last particle
  }

  // Start the outer loop
  for (PartIter1 = StartOuterLoop; PartIter1 != EndOuterLoop; PartIter1++) {

    // If analyzing identical particles, start inner loop at the particle
    // after the current outer loop position, (loops until end)
    if (!partCollection2) {
      StartInnerLoop = PartIter1;
      StartInnerLoop++;
    }

    // If we have two collections - set the first track
    if (mReshuffle != 1 || partCollection2 != nullptr) {
      ThePair->setTrack1(*PartIter1);
    }
    // Start the inner loop
    for (PartIter2 = StartInnerLoop; PartIter2 != EndInnerLoop; PartIter2++) {
      //If we have two collections  - only set the second track
      if (mReshuffle != 1 || partCollection2 != nullptr) {
	ThePair->setTrack2(*PartIter2);
      } else {
	// Swap between first and second particles to avoid biased ordering
	ThePair->setTrack1(swpart ? *PartIter2 : *PartIter1);
	ThePair->setTrack2(swpart ? *PartIter1 : *PartIter2);
	swpart = !swpart;
      }
      // Check if the pair passes the cut
      bool tmpPassPair = mPairCut->pass(ThePair);

      mPairCut->fillCutMonitor(ThePair, tmpPassPair);

      // If pair passes cut, loop over CF's and add pair to real/mixed
      if (tmpPassPair) {

	// Iterate over correlation functions
	for (CorrFctnIter = mCorrFctnCollection->begin();
	     CorrFctnIter != mCorrFctnCollection->end(); CorrFctnIter++) {

	  MpdFemtoBaseCorrFctn* CorrFctn = *CorrFctnIter;
	  if (type == "real") {
	    CorrFctn->addRealPair(ThePair);
	  } else if (type == "mixed") {
	    CorrFctn->addMixedPair(ThePair);
	  } else {
	    std::cout << "Problem with pair type, type = " << type.c_str() << std::endl;
	  }

	} //for (CorrFctnIter=mCorrFctnCollection->begin() ...
      } //if (mPairCut->Pass(ThePair))
    } //for (PartIter2 = StartInnerLoop; PartIter2 != EndInnerLoop; PartIter2++)
  } //for ( PartIter1 = StartOuterLoop; PartIter1 != EndOuterLoop; PartIter1++)

    // We are done with the pair
  delete ThePair;
}

//_________________
void MpdFemtoAnalysis::setReshuffle(unsigned int type) {
  if (type > 2) {
    std::cout << "[WARNING] MpdFemtoAnalsysis::setReshuffle -- Bad reshuffling value \n Resetting...\n";
    mReshuffle = 1;
  } else {
    if (type == 0) {
      mReshuffle = 0;
      return;
    }
    if (type == 1) {
      mReshuffle = 1;
      return;
    }
    if (type == 2) {
      mReshuffle = 2;
      return;
    }
  }
}

//_________________
void MpdFemtoAnalysis::eventBegin(const MpdFemtoEvent* ev) {

  // Perform initialization operations at the beginning of the event processing
  //std::cout << " MpdFemtoAnalysis::EventBegin(const MpdFemtoEvent* ev) " << std::endl;
  mFirstParticleCut->eventBegin(ev);
  mSecondParticleCut->eventBegin(ev);
  mPairCut->eventBegin(ev);
  for (MpdFemtoCorrFctnIterator iter = mCorrFctnCollection->begin();
       iter != mCorrFctnCollection->end(); iter++) {
    (*iter)->eventBegin(ev);
  }
}

//_________________
void MpdFemtoAnalysis::eventEnd(const MpdFemtoEvent* ev) {

  // Finish operations at the end of event processing
  mFirstParticleCut->eventEnd(ev);
  mSecondParticleCut->eventEnd(ev);
  mPairCut->eventEnd(ev);
  for (MpdFemtoCorrFctnIterator iter = mCorrFctnCollection->begin();
       iter != mCorrFctnCollection->end(); iter++) {
    (*iter)->eventEnd(ev);
  }
}

//_________________
void MpdFemtoAnalysis::finish() {
  MpdFemtoCorrFctnIterator iter;
  for (iter = mCorrFctnCollection->begin();
       iter != mCorrFctnCollection->end(); iter++) {
    (*iter)->finish();
  }
}

//_________________
void MpdFemtoAnalysis::addEventProcessed() {
  mNeventsProcessed++;
}

//_________________
TList* MpdFemtoAnalysis::listSettings() {

  // Collect settings list
  const TString setting_prefix = "MpdFemtoAnalysis.";

  TList *tListSettings = new TList();

  TList *event_cut_settings = mEventCut->listSettings();
  if (event_cut_settings != nullptr) {
    TListIter next_event_setting(event_cut_settings);
    while (TObject * obj = next_event_setting()) {
      tListSettings->Add(new TObjString(setting_prefix + obj->GetName()));
    }
    delete event_cut_settings;
  } //if (event_cut_settings != nullptr)

  TList *p1_cut_settings = mFirstParticleCut->listSettings();
  if (p1_cut_settings != nullptr) {
    TListIter next_p1_setting(p1_cut_settings);
    while (TObject * obj = next_p1_setting()) {
      tListSettings->Add(new TObjString(setting_prefix + obj->GetName()));
    }
    delete p1_cut_settings;
  } //if (p1_cut_settings != nullptr)

  if (mSecondParticleCut != mFirstParticleCut) {
    TList *p2_cut_settings = mSecondParticleCut->listSettings();
    if (p2_cut_settings != nullptr) {
      TListIter next_p2_setting(p2_cut_settings);
      while (TObject * obj = next_p2_setting()) {
	tListSettings->Add(new TObjString(setting_prefix + obj->GetName()));
      }
      delete p2_cut_settings;
    }
  } //if (fSecondParticleCut != fFirstParticleCut)

  TList *pair_cut_settings = mPairCut->listSettings();
  if (pair_cut_settings != nullptr) {
    TListIter next_pair_setting(pair_cut_settings);
    while (TObject * obj = next_pair_setting()) {
      tListSettings->Add(new TObjString(setting_prefix + obj->GetName()));
    }
    delete pair_cut_settings;
  } //if (pair_cut_settings != nullptr)

  return tListSettings;
}

//_________________________
TList* MpdFemtoAnalysis::getOutputList() {

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
