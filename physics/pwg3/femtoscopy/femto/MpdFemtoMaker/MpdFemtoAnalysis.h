/**
 * \class MpdFemtoAnalsysis
 * \brief An example of the most basic (concrete) analysis.
 *
 * Most other analyses (e.g. MpdFemtoVertexAnalysis) inherit from this one.
 * Provides basic functionality for the analysis. To properly set up the
 * analysis the following steps should be taken:
 *
 * - create particle cuts and add them via setFirstParticleCut and
 *  setSecondParticleCut. If one analyzes identical particle
 *  correlations, the first particle cut must be also the second
 *  particle cut.
 *
 * - create pair cuts and add them via setPairCut
 *
 * - create one or many correlation functions and add them via
 *  addCorrFctn method.
 *
 * - specify how many events are to be strored in the mixing buffer for
 *  background construction
 *
 * Then, when the analysis is run, for each event, the eventBegin is
 * called before any processing is done, then the ProcessEvent is called
 * which takes care of creating real and mixed pairs and sending them
 * to all the registered correlation functions. At the end of each event,
 * after all pairs are processed, eventEnd is called. After the whole
 * analysis finishes (there is no more events to process) finish() is
 * called.
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoAnalysis_h
#define MpdFemtoAnalysis_h

// MpdFemtoMaker headers
// Base classes
#include "MpdFemtoBaseAnalysis.h"
#include "MpdFemtoBasePairCut.h"
#include "MpdFemtoBaseEventCut.h"
#include "MpdFemtoBaseParticleCut.h"
#include "MpdFemtoBaseCorrFctn.h"

// Infrustructure classes
#include "MpdFemtoCorrFctnCollection.h"
#include "MpdFemtoPicoEventCollection.h"
#include "MpdFemtoParticleCollection.h"
#include "MpdFemtoPicoEvent.h"

// ROOT headers
#include "TList.h"

// Forward declaration
class MpdFemtoPicoEventCollectionVectorHideAway;

//_________________
class MpdFemtoAnalysis : public MpdFemtoBaseAnalysis {
 public:
  /// Construct with default parameters
  ///
  /// All pointer members are initialized to NULL except for the correlation
  /// function collection (mCorrFctnCollection) and the mixing buffer
  /// (mMixingBuffer) which are created with default parameters.
  MpdFemtoAnalysis();

  /// Copy parameters from another analysis.
  ///
  /// All parameters are copied and cuts & correlation functions are cloned.
  /// A new (empty) mixing buffer is created and the number of events processed
  /// (mNeventsProcessed) is set to 0. The EventCollectionHideAway is NOT
  /// copied, and it's up to the subclass to clone if neccessary.
  MpdFemtoAnalysis(const MpdFemtoAnalysis&);
  /// Copy constructor
  MpdFemtoAnalysis& operator=(const MpdFemtoAnalysis&);
  /// Default destructor
  virtual ~MpdFemtoAnalysis();

  //
  // Setters and getters
  //

  /// Return pointer to a pair cut
  virtual MpdFemtoBasePairCut* pairCut();
  /// Return pointer to an event cut
  virtual MpdFemtoBaseEventCut* eventCut();
  /// Return pointer to a first particle cut
  virtual MpdFemtoBaseParticleCut* firstParticleCut();
  /// Return pointer to a first particle cut
  virtual MpdFemtoBaseParticleCut* secondParticleCut();

  /// Return pointer to the correlation function collection
  MpdFemtoCorrFctnCollection* corrFctnCollection();
  
  /// Access to CFs within the collection
  virtual MpdFemtoBaseCorrFctn* corrFctn(int n);
  
  /// Add correlation function to the analysis
  void addCorrFctn(MpdFemtoBaseCorrFctn* cf) {
    mCorrFctnCollection->push_back(cf);
    cf->setAnalysis((MpdFemtoBaseAnalysis*)this);
  }

  /// Set pair cut
  void setPairCut(MpdFemtoBasePairCut* x) {
    mPairCut = x;
    x->setAnalysis((MpdFemtoBaseAnalysis*)this);
  }
  
  /// Set event cut
  void setEventCut(MpdFemtoBaseEventCut* x) {
    mEventCut = x;
    x->setAnalysis((MpdFemtoBaseAnalysis*)this);
  }
  
  /// Set first particle cut
  void setFirstParticleCut(MpdFemtoBaseParticleCut* x) {
    mFirstParticleCut = x;
    x->setAnalysis((MpdFemtoBaseAnalysis*)this);
  }
  
  /// Set second particle cut
  void setSecondParticleCut(MpdFemtoBaseParticleCut* x) {
    mSecondParticleCut = x;
    x->setAnalysis((MpdFemtoBaseAnalysis*)this);
  }

  /// Set minimal size of the particle collection
  void setMinSizePartCollection(unsigned int& minSize) {
    mMinSizePartCollection = minSize;
  }
  
  /// Set Reshuffling mode for particles
  /// \par 0  as is (do not swap of reshuffle particles)
  /// \par 1  swap particles from pair
  /// \par 2  reshuffle using Fisher Yates algorythm
  void setReshuffle(unsigned int type = 1);
  
  /// Set verbose mode
  void setVerboseMode(const bool& isVerbose) {
    mVerbose = isVerbose;
  }


  /// Return size of the event buffer to mix
  unsigned int numEventsToMix() {
    return mNumEventsToMix;
  }
  
  /// Set number of events to mix
  void setNumEventsToMix(const unsigned int& nmix) {
    mNumEventsToMix = nmix;
  }
  
  /// Return pointer ot the current event
  MpdFemtoPicoEvent* currentPicoEvent() {
    return mPicoEvent;
  }
  
  /// Return pointer to the current mixing buffer (particle collection)
  MpdFemtoPicoEventCollection* mixingBuffer() {
    return mMixingBuffer;
  }
  
  /// If mixing buffer is full
  bool mixingBufferFull() {
    return ( mMixingBuffer->size() >= mNumEventsToMix);
  }
  
  /// If first and second particle are identical
  bool analyzeIdenticalParticles() {
    return (mFirstParticleCut == mSecondParticleCut);
  }
  
  ///Returns particles reshuffling type
  char reshuffle() {
    return mReshuffle;
  }
  /// Returns report of reshuffle settings
  MpdFemtoString reshReport();
  /// Returns reports of all cuts applied and correlation functions being done
  virtual MpdFemtoString report();
  /// Return list of cut settings for the analysis
  virtual TList* listSettings();
  /// Return a TList of objects to be written as output
  virtual TList* getOutputList();

  /// Initialization code run at the beginning of processing an event
  ///
  /// This is implemented by calling EventBegin for each member cut
  /// and correlation function
  virtual void eventBegin(const MpdFemtoEvent*);

  /// Bulk of analysis code
  ///
  /// This functions begins by calling EventBegin. If the event passes the
  /// event cut, pairs are made from the particles passing their respective
  /// cuts. The pairs are passed to each correlation function's AddRealPair
  /// method. Pairs made between particles in this event and events in the
  /// mixing buffer, are passed to the correlation functions' AddMixedPair
  /// method. The event is then added to the mixing buffer. The EventEnd() is
  /// called exactly once upon exiting this function.
  virtual void processEvent(const MpdFemtoEvent*);

  /// Cleanup code after processing each event
  ///
  /// Calls EventEnd for each member cut and correlation function.
  virtual void eventEnd(const MpdFemtoEvent*);

  /// Returns number of events which have been passed to processEvent.
  int nEventsProcessed() {
    return mNeventsProcessed;
  }

  /// Finish the analysis
  virtual void finish();

  friend class MpdFemtoLikeSignAnalysis;

 protected:

  /// Increment fNeventsProcessed - is this method neccessary?
  void addEventProcessed();

  /// Build pairs, check pair cuts, and call CFs' AddRealPair() or
  /// AddMixedPair() methods. If no second particle collection is
  /// specfied, make pairs within first particle collection.
  ///
  /// \param type Either the string "real" or "mixed", specifying which method
  ///             to call (AddRealPair or AddMixedPair)
  void makePairs(const char* type, MpdFemtoParticleCollection*, MpdFemtoParticleCollection* p2 = 0);

  /// Mixing Buffer used for Analyses which wrap this one
  MpdFemtoPicoEventCollectionVectorHideAway* mPicoEventCollectionVectorHideAway; //!

  /// Pair cut
  MpdFemtoBasePairCut* mPairCut;
  /// Correlation function collection
  MpdFemtoCorrFctnCollection* mCorrFctnCollection;
  /// Event cut
  MpdFemtoBaseEventCut* mEventCut;
  /// The first particle cut
  MpdFemtoBaseParticleCut* mFirstParticleCut;
  /// The second particle cut
  MpdFemtoBaseParticleCut* mSecondParticleCut;
  /// Event collection (mixing buffer)
  MpdFemtoPicoEventCollection* mMixingBuffer;
  /// Pointer to the event
  MpdFemtoPicoEvent* mPicoEvent; //!

  /// How many "previous" events get mixed with this one, to make background
  unsigned int mNumEventsToMix;
  /// Don't use event if it has fewer than this many particles passing ParticleCuts default 0
  unsigned int mNeventsProcessed;
  /// Minimum # particles in ParticleCollection
  unsigned int mMinSizePartCollection;
  /// Particles Reshuffling
  ///0 - no reshuffling
  ///1 - swapping particles in pairs
  ///2 - Fisher–Yates shuffling algorithm for particle collection

  /// default: mReshuffle=1
  char mReshuffle;
  /// Print info
  bool mVerbose;

  ClassDef(MpdFemtoAnalysis, 0)
};

//_________________
inline MpdFemtoBasePairCut* MpdFemtoAnalysis::pairCut() {
  return mPairCut;
}

//_________________
inline MpdFemtoBaseEventCut* MpdFemtoAnalysis::eventCut() {
  return mEventCut;
}

//_________________
inline MpdFemtoBaseParticleCut* MpdFemtoAnalysis::firstParticleCut() {
  return mFirstParticleCut;
}

//_________________
inline MpdFemtoBaseParticleCut* MpdFemtoAnalysis::secondParticleCut() {
  return mSecondParticleCut;
}

//_________________
inline MpdFemtoCorrFctnCollection* MpdFemtoAnalysis::corrFctnCollection() {
  return mCorrFctnCollection;
}

#endif // #define MpdFemtoAnalysis_h
