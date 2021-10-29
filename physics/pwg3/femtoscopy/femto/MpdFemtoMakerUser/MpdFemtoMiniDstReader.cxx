//
// Reader for the MpdMiniDst format
//

// MpdFemtoMaker headers
#include "MpdFemtoBaseEventCut.h"
#include "MpdFemtoBaseTrackCut.h"
#include "MpdFemtoBaseV0Cut.h"
#include "MpdFemtoBaseKinkCut.h"
#include "MpdFemtoBaseXiCut.h"

// MpdFemtoMakerUser headers
#include "MpdFemtoMiniDstReader.h"

// ROOT headers
#include "TMath.h"
#include "TMatrixDSymEigen.h"
#include "TVectorT.h"
#include "TMatrixTLazy.h"

#ifdef __ROOT__
ClassImp(MpdFemtoMiniDstReader)
#endif

//_________________
MpdFemtoMiniDstReader::MpdFemtoMiniDstReader() :
MpdFemtoBaseEventReader(),
  mMiniDstReader(nullptr),
  mHbtEvent(nullptr) {
  // Default Constructor
  
  // 0 - means good
  mReaderStatus = 0;
}

//_________________
MpdFemtoMiniDstReader::MpdFemtoMiniDstReader(MpdMiniDstReader *miniDstReader, int debug) :
  MpdFemtoBaseEventReader(),
  mMiniDstReader(miniDstReader),
  mHbtEvent(nullptr),
  mEventsPassed(0) {
  // Parametrized constructor
  mDebug = debug;
  // 0 - means good
  mReaderStatus = 0;
}

//_________________
MpdFemtoMiniDstReader::MpdFemtoMiniDstReader(const MpdFemtoMiniDstReader& copy) :
  MpdFemtoBaseEventReader(),
  mMiniDstReader(copy.mMiniDstReader),
  mHbtEvent(nullptr),
  mEventsPassed(0) {
  // Copy constructor

  // 0 - means good
  mReaderStatus = 0;
}

//_________________
MpdFemtoMiniDstReader& MpdFemtoMiniDstReader::operator=(const MpdFemtoMiniDstReader& copy) {
  // Assignment operator
  if (this != &copy) {
    mMiniDstReader = copy.mMiniDstReader;
    mHbtEvent = nullptr;
  } // if (this != &copy)

  // 0 - means good
  mReaderStatus = 0;
  mEventsPassed = 0;

  return *this;
}

//_________________
MpdFemtoMiniDstReader::~MpdFemtoMiniDstReader() {
  if(mMiniDstReader)   { delete mMiniDstReader; mMiniDstReader = nullptr; }
  if(mHbtEvent)        { delete mHbtEvent; mHbtEvent = nullptr; }
}

//_________________
MpdFemtoString MpdFemtoMiniDstReader::report() {
  // Make report
  MpdFemtoString repstr = "\nMpdFemtoString MpdFemtoMiniDstReader::report() - reporting\n";

  repstr += "---> Event cut in the reader: ";
  if (mEventCut) {
    repstr += mEventCut->report();
  }
  else {
    repstr += "\tNONE";
  }

  repstr += "\n---> Track cut in the reader: ";
  if (mTrackCut) {
    repstr += mTrackCut->report();
  }
  else {
    repstr += "\tNONE";
  }

  repstr += "\n---> V0 cut in the reader: ";
  if (mV0Cut) {
    repstr += mV0Cut->report();
  }
  else {
    repstr += "\tNONE";
  }

  repstr += "\n---> V0 cut in the reader: ";
  if (mV0Cut) {
    repstr += mV0Cut->report();
  }
  else {
    repstr += "\tNONE";
  }

  repstr += "\n---> Kink cut in the reader: ";
  if (mKinkCut) {
    repstr += mKinkCut->report();
  }
  else {
    repstr += "\tNONE";
  }

  repstr += "\n---> Xi cut in the reader: ";
  if (mXiCut) {
    repstr += mXiCut->report();
  }
  else {
    repstr += "\tNONE";
  }

  repstr += "\n";
  return repstr;
}

//_________________
MpdFemtoEvent* MpdFemtoMiniDstReader::returnHbtEvent() {
  // Read miniEvent, create MpdFemtoEvent and return it to the MpdFemtoAnalysis

  // Clean-up mHbtEvent from previous reading
  mHbtEvent = nullptr;

  // Check that MpdMiniDstReader exists
  if (!mMiniDstReader) {
    std::cout << "[ERROR] MpdFemtoEvent* MpdFemtoMiniDstReader::returnHbtEvent() - no MpdMiniDstReader is provided" << std::endl;
    mReaderStatus = 1;
    return mHbtEvent;
  }

  // Read (load) event
  Bool_t isGood = mMiniDstReader->readMiniEvent(mEventsPassed);
  // Increment counter
  mEventsPassed++;

  // Check that miniDst exists
  MpdMiniDst* miniDst = (MpdMiniDst*)mMiniDstReader->miniDst();
  if (!miniDst) {
    std::cout << "[ERROR] MpdFemtoEvent* MpdFemtoMiniDstReader::returnHbtEvent() - no miniDst is in MpdMiniDstReader" << std::endl;
    mReaderStatus = 1;
    return mHbtEvent;
  }

  // Check that MpdMiniEvent exists
  if ( !miniDst->event() || !miniDst->numberOfTracks() ) {// CAN'T UNDERSTAND WHY THERE IS NO "!" BEFORE SECOND CONDITION or <=0
    std::cout << "[WARNING] MpdFemtoEvent* MpdFemtoMiniDstReader::returnHbtEvent()"
              << " - MpdMiniEvent does not exist or particles missing"
              << std::endl;
    mReaderStatus = 1;
  }
  else {

    // Create new MpdFemtoEvent
    mHbtEvent = new MpdFemtoEvent();

    // Retrieve MpdMiniEvent
    MpdMiniEvent *miniEvent = (MpdMiniEvent*)miniDst->event();

    float bField = miniEvent->bField();

    mHbtEvent->setEventId( miniEvent->eventId() );
    mHbtEvent->setRunNumber( miniEvent->runId() );
    mHbtEvent->setBField( bField );
    mHbtEvent->setRefMult( miniEvent->refMult() );
    mHbtEvent->setRefMultPos( miniEvent->refMultPos() );
    mHbtEvent->setGRefMult( miniEvent->grefMult() );

    mHbtEvent->setCent16( -1 );
    mHbtEvent->setRefMultCorr(-1);
    mHbtEvent->setRefMultCorrWeight( 1 );

    mHbtEvent->setBTofTrayMult( miniEvent->btofTrayMultiplicity() );
    mHbtEvent->setNumberOfBTofMatched( miniEvent->nBTOFMatch() );
    mHbtEvent->setNumberOfBEMCMatched( miniEvent->nBECalMatch() );
    // Could be recalculated, but skipped for now
    //mHbtEvent->setZdcAdcEast( miniEvent->ZdcSumAdcEast() );
    //mHbtEvent->setZdcAdcWest( miniEvent->ZdcSumAdcWest() );

    mHbtEvent->setZdcCoincidenceRate( miniEvent->fhcalX() );

    mHbtEvent->setPrimaryVertex( miniEvent->primaryVertex() );
    mHbtEvent->setVpdVz( miniEvent->vzFfd() );
    mHbtEvent->setTriggerIds( miniEvent->triggerIds() );

    // Loop over miniTracks in the event
    for (unsigned int iTrk=0; iTrk<miniDst->numberOfTracks(); iTrk++) {

      // Retrieve i-th miniTrack
      MpdMiniTrack *miniTrack = (MpdMiniTrack*)miniDst->track(iTrk);
      if (!miniTrack) continue;

      // Create new instance of the MpdFemtoTrack
      MpdFemtoTrack *hbtTrack = new MpdFemtoTrack();
      hbtTrack->setId( miniTrack->id() );
      hbtTrack->setNHits( miniTrack->nHits() * miniTrack->charge() );
      //hbtTrack->setNHitsPossible( miniTrack->nHitsPoss() );
      //hbtTrack->setNHitsDedx( miniTrack->nHitsDedx() );
      hbtTrack->setChi2( miniTrack->chi2() );
      hbtTrack->setDedxFromKeV( miniTrack->dEdx() );
      hbtTrack->setNSigmaElectron( miniTrack->nSigmaElectron() );
      hbtTrack->setNSigmaPion( miniTrack->nSigmaPion() );
      hbtTrack->setNSigmaKaon( miniTrack->nSigmaKaon() );
      hbtTrack->setNSigmaProton( miniTrack->nSigmaProton() );
      hbtTrack->setPidProbElectron( 0.5 );
      hbtTrack->setPidProbPion( 0.5 );
      hbtTrack->setPidProbKaon( 0.5 );
      hbtTrack->setPidProbProton( 0.5 );
      hbtTrack->setDca( miniTrack->origin().X() - miniEvent->primaryVertex().X(),
                        miniTrack->origin().Y() - miniEvent->primaryVertex().Y(),
                        miniTrack->origin().Z() - miniEvent->primaryVertex().Z() );
      if ( miniTrack->isPrimary() ) {
        hbtTrack->setP( miniTrack->pMom() );
      }
      else {
        hbtTrack->setP( 0., 0., 0. );
      }
      hbtTrack->setGlobalP( miniTrack->gMom() );
      hbtTrack->setPrimaryVertex( miniEvent->primaryVertex() );
      hbtTrack->setMagneticField( bField );
      hbtTrack->setTopologyMap( miniTrack->hitMap() );
      if ( miniTrack->isBTofTrack() ) {
        MpdMiniBTofPidTraits *trait = miniDst->btofPidTraits( miniTrack->bTofPidTraitsIndex() );
        if (!trait) {
          std::cout << "[WARNING] MpdFemtoEvent* MpdFemtoMiniDstReader::returnHbtEvent()"
                    << " -  No MpdMiniBTofPidTraits was found by the index" << std::endl;
          hbtTrack->setBeta( -999. );
        }
        else {
          hbtTrack->setBeta( trait->beta() );
        }
      }
      else {
        hbtTrack->setBeta( -999. );
      }

      // Hidden info should be used for the MC data
      MpdFemtoModelHiddenInfo *hiddenInfo = new MpdFemtoModelHiddenInfo();
      
      Int_t mcTrackId = -1; // filled later with real index of corresp. Monte-Carlo track ...
      
      if (hbtTrack->IsPrimary()) {        
         mcTrackId = miniTrack->mcTrackIndex();
      
        if (mcTrackId != -1) {
         MpdMiniMcTrack *mcTrack = miniDst->mcTrack(mcTrackId);
         
         // Setting momentum ...
         hiddenInfo->setTrueMomentum(mcTrack->px(), mcTrack->py(), mcTrack->pz());
         
         // Setting emission point ...
         hiddenInfo->setEmissionPoint(mcTrack->x(), mcTrack->y(), mcTrack->z(), mcTrack->t());
         
         // Setting pdg ..
         hiddenInfo->setPdgPid(mcTrack->pdgId());
         hiddenInfo->setOrigin(0);
        }
      }

      if (hiddenInfo)
        hbtTrack->setHiddenInfo( hiddenInfo );

      // Check if a front-loaded cut exists. The mTrackCut
      // is inherited from MpdFemtoBaseReader
      if ( mTrackCut && !mTrackCut->pass( hbtTrack ) ) {
        delete hbtTrack;
        continue;
      } // if ( mTrackCut )

      mHbtEvent->trackCollection()->push_back( hbtTrack );

    } // for (int iTrk=0; iTrk<miniDst->numberOfTracks(); iTrk++)

    // Check front-loaded cuts here, because now all event information
    // should be filled making event cut on it possible. The mEventCut
    // is inherited from MpdFemtoBaseReader
    if (mEventCut && !mEventCut->pass(mHbtEvent) ) {
      delete mHbtEvent;
      mHbtEvent = nullptr;
      // return mHbtEvent; // Will be deleted 6 lines below
    } // if (mEventCut)

  } // else

  return mHbtEvent;
}
