// MpdRoot headers
#include "MpdEvent.h"
#include "MpdVertex.h"
#include "MpdTrack.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdKalmanFilter.h"
#include "MpdKalmanGeoScheme.h"
#include "MpdTofHit.h"
#include "MpdEmcDigitKI.h"
#include "MpdEmcClusterKI.h"
#include "MpdTofMatchingData.h"
#include "MpdZdcDigi.h"
#include "MpdMCTrack.h"
#include "MpdGenTrack.h"

// MiniDst headers
#include "MpdMiniDstFillTask.h"
#include "MpdMiniDst.h"
#include "MpdMiniEvent.h"
#include "MpdMiniTrack.h"
#include "MpdMiniBTofHit.h"
#include "MpdMiniBTofPidTraits.h"
#include "MpdMiniBECalCluster.h"
#include "MpdMiniTrackCovMatrix.h"
#include "MpdMiniFHCalHit.h"
#include "MpdMiniMcEvent.h"
#include "MpdMiniMcTrack.h"
#include "MpdMiniMessMgr.h"
#include "MpdMiniArrays.cxx"

//_________________

MpdMiniDstFillTask::MpdMiniDstFillTask()
   : fIsUseCovMatrix(kTRUE), fIsUseECal(kTRUE), fStoreMatched(kFALSE), fStoreParents(kFALSE), fTrivDecay(kTRUE),
     fEventHeaders(nullptr), fEvents(nullptr), fVertices(nullptr), fTpcTracks(nullptr), fTofHits(nullptr),
     fTofMatching(nullptr), fMCTracks(nullptr), fGenTracks(nullptr), fEmcClusters(nullptr), fZdcDigits(nullptr),
     fMiniDst(nullptr), fBField(0.5), // mag. field in T (MUST be changed to the real magnetic field
     fOutputFile(nullptr), fTTree(nullptr), fSplit(99), fCompression(9), fBufferSize(65536 * 4), fMiniArrays(nullptr),
     fNSigmaDedxEstimator(0)
{
   // Default constructor
   /* empty */
}

//_________________

MpdMiniDstFillTask::MpdMiniDstFillTask(TString name) : MpdMiniDstFillTask()
{
   fMiniDst = new MpdMiniDst();
   // Standard constructor
   TString truncatedInFileName = name;
   Ssiz_t  lastOccurence       = name.Last('/');
   if (lastOccurence != kNPOS) {
      truncatedInFileName = name((lastOccurence + 1), name.Length());
   }
   fOutputFileName = truncatedInFileName.ReplaceAll(".root", ".MiniDst.root");

   if (!fMcTrk2MiniMcTrk.empty()) fMcTrk2MiniMcTrk.clear();
   if (!fMcTrk2EcalCluster.empty()) fMcTrk2EcalCluster.clear();

   streamerOff();
   createArrays();
}

//_________________

MpdMiniDstFillTask::~MpdMiniDstFillTask()
{
   // Destructor
   if (!fMcTrk2MiniMcTrk.empty()) fMcTrk2MiniMcTrk.clear();
   if (!fMcTrk2EcalCluster.empty()) fMcTrk2EcalCluster.clear();
   if (fMiniDst) delete fMiniDst;
}

//_________________

InitStatus MpdMiniDstFillTask::Init()
{
   // Output name must exist
   if (fOutputFileName.IsNull()) {
      return kERROR;
   }

   // Creating output tree with miniDST
   fOutputFile = new TFile(fOutputFileName.Data(), "RECREATE");
   // Inform about the creation
   LOG_INFO << " Output file: " << fOutputFileName.Data() << " created." << endm;
   // Set compression level
   fOutputFile->SetCompressionLevel(fCompression);
   int bufsize = fBufferSize;
   if (fSplit) {
      bufsize /= 4;
   }
   // Create TTree
   fTTree = new TTree("MiniDst", "MpdMiniDst", fSplit);
   fTTree->SetAutoSave(1000000);

   // Create arrays
   for (Int_t i = 0; i < MpdMiniArrays::NAllMiniArrays; ++i) {
      fTTree->Branch(MpdMiniArrays::miniArrayNames[i], &fMiniArrays[i], bufsize, fSplit);
   }

   // Initialize FairRoot manager
   FairRootManager *ioman = FairRootManager::Instance();
   if (!ioman) {
      std::cout << "[ERROR] MpdMiniDstFillTask::Init - Not FairRootManager instance has been found" << std::endl;
      exit(0);
   }

   // Reading all possible input branches ...

   // Get MC event header
   fEventHeaders = (FairMCEventHeader *)ioman->GetObject("MCEventHeader.");
   if (!fEventHeaders) {
      std::cout << "[WARNING] MpdMiniDstFillTask::Init - No MCEventHeader has been found" << std::endl;
   }

   // Get information on event
   fEvents = (MpdEvent *)ioman->GetObject("MPDEvent.");
   if (!fEvents) {
      std::cout << "[WARNING] MpdMiniDstFillTask::Init - No MpdEvent has been found" << std::endl;
   }

   // Get information on reconstructed vertex
   fVertices = (TClonesArray *)ioman->GetObject("Vertex");
   if (!fVertices) {
      std::cout << "[WARNING] MpdMiniDstFillTask::Init - No Vertex has been found" << std::endl;
   }

   // Retrieve information about tracks in TPC
   fTpcTracks = (TClonesArray *)ioman->GetObject("TpcKalmanTrack");
   if (!fTpcTracks) {
      std::cout << "[WARNING] MpdMiniDstFillTask::Init - No TpcKalmanTrack has been found" << std::endl;
   }

   // Retrieve information about TOF hits
   fTofHits = (TClonesArray *)ioman->GetObject("TOFHit");
   if (!fTofHits) {
      std::cout << "[WARNING] MpdMiniDstFillTask::Init - No TOFHit has been found" << std::endl;
   }

   // Get information about TPC tracks that match TOF
   fTofMatching = (TClonesArray *)ioman->GetObject("TOFMatching");
   if (!fTofMatching) {
      std::cout << "[WARNING] MpdMiniDstFillTask::Init - No TOFMatching has been found" << std::endl;
   }

   // Retrieve MC Tracks
   fMCTracks = (TClonesArray *)ioman->GetObject("MCTrack");
   if (!fMCTracks) {
      std::cout << "[WARNING] MpdMiniDstFillTask::Init - No MCTracks have been found" << std::endl;
   }

   // Retrieve generator (primary) tracks
   fGenTracks = (TClonesArray *)ioman->GetObject("GenTracks");
   if (!fGenTracks) {
      std::cout << "[WARNING] MpdMiniDstFillTask::Init - No GenTracks have been found" << std::endl;
   }

   // Retrieve information about barrel ECal digits
   // fEmcDigits = (TClonesArray*) ioman->GetObject("EmcDigit");
   // if (!fEmcDigits) {
   //   std::cout << "[WARNING] MpdMiniDstFillTask::Init - No EmcDigits have been found"
   // 	      << std::endl;
   // }

   if (fIsUseECal) {
      // Retrieve information about barrel ECal clusters
      fEmcClusters = (TClonesArray *)ioman->GetObject("EmcCluster");
      if (!fEmcClusters) {
         std::cout << "[WARNING] MpdMiniDstFillTask::Init - No EmcClusters have been found" << std::endl;
      }
   } // if ( fIsUseECal )

   // Retrieve information about FHCal
   fZdcDigits = (TClonesArray *)ioman->GetObject("ZdcDigi");
   if (!fZdcDigits) {
      std::cout << "[WARNING] MpdMiniDstFillTask::Init - No ZdcDigits have been found" << std::endl;
   }

   return kSUCCESS;
}

//_________________

void MpdMiniDstFillTask::Exec(Option_t *option)
{

   // Clear container at the beginning of each event
   if (!fMcTrk2MiniMcTrk.empty()) fMcTrk2MiniMcTrk.clear();
   if (!fMcTrk2EcalCluster.empty()) fMcTrk2EcalCluster.clear();

   Bool_t isGood = isGoodEvent();

   // Convert only events that are normally reconstructed
   if (isGood) {

      fillEvent();
      if (fZdcDigits) fillFHCalHits();
      fillBTofHits();
      findMatchedMcTracks();
      fillMcTracks();
      if (fIsUseECal) fillECalClusters();
      fillTracks();

      if (fVerbose != 0) {
         fMiniDst->printTracks();
      }

      // Fill TTree
      fTTree->Fill();
   } // if (isGood)
   else {
      std::cout << "[WARNING] MpdMiniDstFillTask::Exec - Bad event reconstuction. Skipping event" << std::endl;
   }
}

//_________________

Bool_t MpdMiniDstFillTask::isGoodEvent()
{

   // Bool_t isGoodAmount = kFALSE;
   // Bool_t isGoodOrder = kTRUE;

   TClonesArray *glTracks      = (TClonesArray *)fEvents->GetGlobalTracks();
   Int_t         nGlobalTracks = glTracks->GetEntriesFast();
   Int_t         nKalmanTracks = fTpcTracks->GetEntriesFast();

   if (nGlobalTracks != nKalmanTracks) return kFALSE;
   // isGoodAmount = kTRUE;

   // Loop over kalman (aka TPC) tracks
   for (Int_t iKTrack = 0; iKTrack < nKalmanTracks; iKTrack++) {
      // Retrieve kalman track
      MpdTpcKalmanTrack *kTrack = (MpdTpcKalmanTrack *)fTpcTracks->UncheckedAt(iKTrack);
      if (!kTrack) continue;
      // Retrieve global track
      MpdTrack *gTrack = (MpdTrack *)glTracks->UncheckedAt(iKTrack);
      if (!glTracks || kTrack->GetTrackID() != gTrack->GetID()) continue;
   } // for (Int_t iKTracks=0; iKTracks<nKalmanTracks; iKTracks++)

   return kTRUE;
}

//_________________

void MpdMiniDstFillTask::Finish()
{
   // Write data to the output file and close it
   if (fOutputFile) {
      fOutputFile->Write();
      fOutputFile->Close();
   }
}

//_________________

void MpdMiniDstFillTask::streamerOff()
{
   MpdMiniEvent::Class()->IgnoreTObjectStreamer();
   MpdMiniTrack::Class()->IgnoreTObjectStreamer();
   MpdMiniBTofHit::Class()->IgnoreTObjectStreamer();
   MpdMiniBECalCluster::Class()->IgnoreTObjectStreamer();
   MpdMiniBTofPidTraits::Class()->IgnoreTObjectStreamer();
   MpdMiniTrackCovMatrix::Class()->IgnoreTObjectStreamer();
   MpdMiniFHCalHit::Class()->IgnoreTObjectStreamer();
   MpdMiniMcEvent::Class()->IgnoreTObjectStreamer();
   MpdMiniMcTrack::Class()->IgnoreTObjectStreamer();
}

//_________________

void MpdMiniDstFillTask::createArrays()
{
   // Create MiniDst arrays
   fMiniArrays = new TClonesArray *[MpdMiniArrays::NAllMiniArrays];

   for (Int_t iArr = 0; iArr < MpdMiniArrays::NAllMiniArrays; iArr++) {
      fMiniArrays[iArr] = new TClonesArray(MpdMiniArrays::miniArrayTypes[iArr], MpdMiniArrays::miniArraySizes[iArr]);
   }

   // Set pointers to the arrays
   fMiniDst->set(fMiniArrays);
}

//_________________

void MpdMiniDstFillTask::fillEvent()
{

   // Fill event information
   TClonesArray *miniEventHeaders = fMiniArrays[MpdMiniArrays::Event];
   miniEventHeaders->Delete();

   // Create new event information
   MpdMiniEvent *miniEvent = new ((*miniEventHeaders)[miniEventHeaders->GetEntriesFast()]) MpdMiniEvent();
   if (!miniEvent) {
      std::cout << "[ERROR] MpdMiniDstFillTask::fillEvent - No MiniEvent has been created" << std::endl;
      // In this case we need to quit (probably completely)
      return;
   }

   // Retrieve information about reconstructed primary vertex
   MpdVertex *vertex = (MpdVertex *)fVertices->UncheckedAt(0);

   // Primary vertex position
   if (vertex) {
      miniEvent->setPrimaryVertexPosition(vertex->GetX(), vertex->GetY(), vertex->GetZ());
   } else {
      // In case that vertex has not been reconstructed
      miniEvent->setPrimaryVertexPosition(-999., -999., -999.);
   }

   // Number of global tracks in the current event
   if (!fEvents) {
      std::cout << "[ERROR] MpdMiniDstFillTask::fillEvent - No event information has been found" << std::endl;
      return;
   }

   // Retrieve global tracks
   TClonesArray *glTracks = (TClonesArray *)fEvents->GetGlobalTracks();
   miniEvent->setNumberOfGlobalTracks(glTracks->GetEntriesFast());

   // Fill McEvent info
   TClonesArray *miniMcEventHeaders = fMiniArrays[MpdMiniArrays::McEvent];
   miniMcEventHeaders->Delete();

   // Retrieve number of primary tracks
   Int_t    nPart = -1;
   Int_t    nColl = -1;
   Double_t rp    = -1;
   // Impact parameter
   Double_t b = fEventHeaders->GetB();
   // Run number
   UInt_t runId = fEventHeaders->GetRunID();
   // Event number
   UInt_t eventId = fEventHeaders->GetEventID();
   // Time
   Double_t time = fEventHeaders->GetT();
   // Fill MC vertex
   TVector3 vtx;
   fEventHeaders->GetVertex(vtx);

   // Create and fill MiniMcEvent
   MpdMiniMcEvent *mcEvent = new ((*miniMcEventHeaders)[miniMcEventHeaders->GetEntriesFast()])
      MpdMiniMcEvent(runId, eventId, rp, b, nPart, nColl, vtx, time);
   mcEvent->setReactionPlaneAngle(fEventHeaders->GetRotZ());
}

//_________________

void MpdMiniDstFillTask::fillFHCalHits()
{

   TClonesArray *miniFHCal = fMiniArrays[MpdMiniArrays::FHCalHit];
   miniFHCal->Delete();

   // Retrieve number of hits (aka channels = detId(1,2)*modId(1-45)*ch(1-42) )
   Int_t nFHCalHits = fZdcDigits->GetEntries();

   // Total number of FHCal modules (aka towers)
   Int_t nFHCalModules = 45;
   // Energy deposition in the modules (towers)
   Float_t eLoss[90] = {};

   Int_t det = 0;
   Int_t mod = 0;

   // Loop over hits
   for (Int_t iHit = 0; iHit < nFHCalHits; iHit++) {

      // Retrieve i-th tower information
      MpdZdcDigi *zdcHit = (MpdZdcDigi *)fZdcDigits->At(iHit);

      if (!zdcHit) continue;

      // FHCal miniHit
      MpdMiniFHCalHit *hit = new ((*miniFHCal)[miniFHCal->GetEntriesFast()]) MpdMiniFHCalHit();

      det = zdcHit->GetDetectorID();
      mod = zdcHit->GetModuleID();

      // Set hit information
      hit->setId(det, mod, zdcHit->GetChannelID());
      hit->setEDep(zdcHit->GetELoss());
      // hit->setAdc( zdcHit->ADC( zdcHit->GetELoss() ) );

      // Measure energy deposition for each module
      eLoss[(det - 1) * nFHCalModules + mod - 1] += zdcHit->GetELoss();
   } // for (Int_t iTower=0; iTower<nFHCalTower; iTower++)

   MpdMiniEvent *miniEvent = (MpdMiniEvent *)fMiniArrays[MpdMiniArrays::Event]->At(0);

   // Fill energy deposition in each tower
   for (Int_t iTow = 0; iTow < (2 * nFHCalModules); iTow++) {
      miniEvent->setFHCalEnergyDepositInModule(iTow, eLoss[iTow]);
   }
}

//_________________

void MpdMiniDstFillTask::fillBTofHits()
{

   // Instantiate MpdMiniBTofHit array
   TClonesArray *miniToF = fMiniArrays[MpdMiniArrays::BTofHit];
   miniToF->Delete();

   if (!fTofHits) return;

   // Loop over TOF hits
   for (Int_t iHit = 0; iHit < fTofHits->GetEntriesFast(); iHit++) {

      // Retrieve TOF hit
      MpdTofHit *tofHit = (MpdTofHit *)fTofHits->UncheckedAt(iHit);

      if (!tofHit) continue;

      MpdMiniBTofHit *miniTofHit = new ((*miniToF)[miniToF->GetEntriesFast()]) MpdMiniBTofHit();

      miniTofHit->setId(tofHit->GetDetectorID());
      miniTofHit->setHitPositionXYZ(tofHit->GetX(), tofHit->GetY(), tofHit->GetZ());
      miniTofHit->setTime(tofHit->GetTime());
   } // for (Int_t iHit = 0; iHit < fTofHits->GetEntriesFast(); iHit++)
}

//_________________

void MpdMiniDstFillTask::fillMcTracks()
{

   //
   // Prepare mapping of MC tracks to reconstructed ones
   // and those used in barrel ECal clusters
   //

   TClonesArray *miniMcTracks = fMiniArrays[MpdMiniArrays::McTrack];
   miniMcTracks->Delete();

   if (!fMCTracks) return;

   // Loop over MC tracks
   for (Int_t iMcTrk = 0; iMcTrk < fMCTracks->GetEntriesFast(); iMcTrk++) {

      // Retrieve MCTrack
      MpdMCTrack *mcTrack = (MpdMCTrack *)fMCTracks->UncheckedAt(iMcTrk);

      // MC track must exist
      if (!mcTrack) continue;

      // Check if MC track is a generator level track
      Bool_t isGenLevelTrack = kFALSE;
      if (mcTrack->GetMotherId() == -1) isGenLevelTrack = kTRUE;

      // Check if generator level or ECal track or need to be stored
      if (!fMapRec[iMcTrk]) continue;

      // Create new MiniMcTrack
      MpdMiniMcTrack *miniMcTrack = new ((*miniMcTracks)[miniMcTracks->GetEntriesFast()]) MpdMiniMcTrack();

      if (!miniMcTrack) {
         std::cout << "[WARNING] MpdMiniDstFillTask::fillTracks - No miniMcTrack has been found" << std::endl;
         continue;
      }

      // Set McTrack information
      miniMcTrack->setId(iMcTrk);
      miniMcTrack->setPdgId(mcTrack->GetPdgCode());
      miniMcTrack->setPx(mcTrack->GetPx());
      miniMcTrack->setPy(mcTrack->GetPy());
      miniMcTrack->setPz(mcTrack->GetPz());
      miniMcTrack->setEnergy(mcTrack->GetEnergy());
      miniMcTrack->setMotherId(mcTrack->GetMotherId());

      // Assume that there is no branch with GenTracks ...
      miniMcTrack->setX(-1.);
      miniMcTrack->setY(-1.);
      miniMcTrack->setZ(-1.);
      miniMcTrack->setT(-1.);

      if (fGenTracks && isGenLevelTrack) {
         for (Int_t iGenTrack = 0; iGenTrack < fGenTracks->GetEntriesFast(); iGenTrack++) {
            MpdGenTrack *genTrack = (MpdGenTrack *)fGenTracks->UncheckedAt(iGenTrack);
            if (genTrack->GetIsUsed()) continue;

            Double_t absMomDiff = TMath::Abs(genTrack->GetMomentum().Mag() - mcTrack->GetP());
            if (absMomDiff < DBL_EPSILON) {
               genTrack->SetIsUsed(kTRUE);

               TLorentzVector spaceTime = genTrack->GetCoordinates();
               miniMcTrack->setX(spaceTime.X());
               miniMcTrack->setY(spaceTime.Y());
               miniMcTrack->setZ(spaceTime.Z());
               miniMcTrack->setT(spaceTime.T());
            } // if (absMomDiff < DBL_EPSILON)
         }    // for (Int_t iGenTrack = 0; iGenTrack < fGenTracks->GetEntriesFast(); iGenTrack++)
      }       // if (fGenTracks)

      // Store indices
      fMcTrk2MiniMcTrk.push_back(std::make_pair(iMcTrk, (UShort_t)miniMcTracks->GetEntriesFast() - 1));
   } // for (Int_t iMcTrk = 0; iMcTrk < fMCTracks->GetEntriesFast(); iMcTrk++)

   // fix mother indexes
   for (int iMiniTrack = 0; iMiniTrack < miniMcTracks->GetEntriesFast(); iMiniTrack++) {

      MpdMiniMcTrack *track = (MpdMiniMcTrack *)miniMcTracks->UncheckedAt(iMiniTrack);
      if (track->getMotherId() < 0) continue; // this is particle without parents
      Int_t newIndex = miniMcIdxFromMcIdx(track->getMotherId());
      if (newIndex >= 0)
         track->setMotherId(newIndex);
      else
         track->setMotherId(-2); // set -2 for particles that are not primary but doesn't have parents
   }
}

//_________________

void MpdMiniDstFillTask::fillECalClusters()
{

   //
   // Fill barrel ECal cluster information. Need to do it here
   // because of the mapping
   //

   // Instantiate MpdMiniBECalCluster array
   TClonesArray *miniEmcClusters = fMiniArrays[MpdMiniArrays::BECalCluster];
   miniEmcClusters->Delete();

   // Check if barrel ECal clusters exist
   if (fEmcClusters) {

      // Variables for smaller and larger dispertion axes
      Float_t lambda1 = 0;
      Float_t lambda2 = 0;

      // Loop over barrel ECal clusters
      for (Int_t iCluster = 0; iCluster < fEmcClusters->GetEntriesFast(); iCluster++) {

         // Retrieve barrel ECal cluster
         MpdEmcClusterKI *cluster = (MpdEmcClusterKI *)fEmcClusters->UncheckedAt(iCluster);

         // Cluster must exist
         if (!cluster) continue;

         // Clear
         lambda1 = 0;
         lambda2 = 0;

         // Create miniCluster
         MpdMiniBECalCluster *miniCluster =
            new ((*miniEmcClusters)[miniEmcClusters->GetEntriesFast()]) MpdMiniBECalCluster();

         // Fill miniCluster info
         miniCluster->setEnergy(cluster->GetE());
         miniCluster->setECore(cluster->GetEcore());
         miniCluster->setECore1p(cluster->GetEcore_1p());
         miniCluster->setECore2p(cluster->GetEcore_2p());
         miniCluster->setTime(cluster->GetTime());
         miniCluster->setXYZ(cluster->GetX(), cluster->GetY(), cluster->GetZ());
         miniCluster->setDPhi(cluster->GetDPhi());
         miniCluster->setDz(cluster->GetDZ());
         miniCluster->setTrackId(miniMcIdxFromMcIdx(cluster->GetTrackIndex()));

         cluster->GetLambdas(lambda1, lambda2);
         miniCluster->setLambdas(lambda1, lambda2);
         miniCluster->setChi2(cluster->GetChi2());
         miniCluster->setNLM(cluster->GetNLM());

         // Fill digit info
         Int_t   id   = -1;
         Float_t eDep = -1.;
         for (Int_t iDigi = 0; iDigi < cluster->GetMultiplicity(); iDigi++) {
            id   = -1;
            eDep = -1.;
            cluster->GetDigitParams(iDigi, id, eDep);
            miniCluster->addDigit(id, eDep);
         } // for (Int_t iDigi = 0; iDigi < cluster->GetMultiplicity(); iDigi++)

         // Fill MC track info
         for (Int_t iTrk = 0; iTrk < cluster->GetNumberOfTracks(); iTrk++) {
            id   = -1;
            eDep = -1.;
            cluster->GetMCTrack(iTrk, id, eDep);
            miniCluster->addMcTrack(miniMcIdxFromMcIdx(id), eDep);
         } // for (Int_t iTrk=0; iTrk < cluster->GetNumberOfTracks(); iTrk++)
      }    // for (Int_t iCluster = 0; iCluster < fEmcClusters->GetEntriesFast(); iCluster++)
   }       // if (fEmcClusters)
}

//_________________

void MpdMiniDstFillTask::fillTracks()
{

   //
   // Fill miniTrack, corresponding covariance matrix
   // and BTof matching information
   //

   // Reconstructed tracks
   TClonesArray *miniTracks = fMiniArrays[MpdMiniArrays::Track];
   miniTracks->Delete();
   // Convariance matrices of reconstructed tracks
   TClonesArray *miniTrackCovMatrices = fMiniArrays[MpdMiniArrays::TrackCovMatrix];
   miniTrackCovMatrices->Delete();
   // Create TOF-matching information
   TClonesArray *miniBTofTraits = fMiniArrays[MpdMiniArrays::BTofPidTraits];
   miniBTofTraits->Delete();

   // Retrieve global tracks from event
   TClonesArray *glTracks = (TClonesArray *)fEvents->GetGlobalTracks();

   // Reconstructed primary vertex in event
   MpdVertex *vtx = (MpdVertex *)fVertices->UncheckedAt(0);

   // Retrieve primary track indices in the fTpcTracks array that were used
   // for the primary vertex reconstruction
   TArrayI *ind = vtx->GetIndices();

   std::vector<Int_t> indices;
   for (Int_t iEle = 0; iEle < ind->GetSize(); iEle++) {
      indices.push_back(ind->At(iEle));
   }

   // Reference multiplicities
   Int_t                 grefMult           = 0;
   Int_t                 refMultPos         = 0;
   Int_t                 refMultNeg         = 0;
   Int_t                 refMultHalfPosEast = 0;
   Int_t                 refMultHalfNegEast = 0;
   Int_t                 refMultHalfPosWest = 0;
   Int_t                 refMultHalfNegWest = 0;
   Int_t                 refMult2PosEast    = 0;
   Int_t                 refMult2NegEast    = 0;
   Int_t                 refMult2PosWest    = 0;
   Int_t                 refMult2NegWest    = 0;
   Int_t                 nTofMatched        = 0;
   Int_t                 nBECalMatched      = 0;
   std::vector<Double_t> nSigma;

   // Fill global and primary track information. All tracks originate
   // from kalman track, so loop over them and copy information
   for (Int_t iKalmanTrk = 0; iKalmanTrk < fTpcTracks->GetEntriesFast(); iKalmanTrk++) {

      // Retrieve i-th TpcKalmanTrack
      MpdTpcKalmanTrack *kalTrack = (MpdTpcKalmanTrack *)fTpcTracks->UncheckedAt(iKalmanTrk);

      // Skip non-existing tracks
      if (!kalTrack) continue;

      // Create miniTrack
      MpdMiniTrack *miniTrack = new ((*miniTracks)[miniTracks->GetEntriesFast()]) MpdMiniTrack();

      // Clean vector for each track
      if (!nSigma.empty()) nSigma.clear();

      // Set track parameters
      miniTrack->setId(kalTrack->GetTrackID());
      miniTrack->setChi2(kalTrack->GetChi2());
      miniTrack->setNHits(kalTrack->GetNofHits() * kalTrack->Charge());

      // Fill track covariance matrix if needed
      if (fIsUseCovMatrix) {
         // Create miniTrack covariance matrix
         MpdMiniTrackCovMatrix *miniTrackCovMatrix =
            new ((*miniTrackCovMatrices)[miniTrackCovMatrices->GetEntriesFast()]) MpdMiniTrackCovMatrix();
         fillCovMatrix(kalTrack, miniTrackCovMatrix);
      }

      // Let's find primary tracks from the whole set of reconstructed
      // tracks from the current event
      Bool_t isPrimary = kFALSE;

      for (auto it : indices) {
         if (TMath::Abs(it - iKalmanTrk) == 0) {
            isPrimary = kTRUE;
            break;
         }
      }

      // If track is primary then fill primary momentum
      if (isPrimary) {

         // Refit track to primary vertex
         refit2Vp(miniTrack, iKalmanTrk, vtx);

         //
         // TODO: Here we should check TOF-matching of the primary track
         //

         // Estimate reference multiplicities
         if (miniTrack->nHits() > 15) {

            TVector3 pMom   = miniTrack->pMom();
            Int_t    charge = (miniTrack->charge() > 0) ? 1 : -1;

            // Central region
            if (TMath::Abs(pMom.Eta()) < 0.5) {
               (charge > 0) ? refMultPos++ : refMultNeg++;
            }

            // Halfs of the TPC
            if (pMom.Eta() > 0) {
               (charge > 0) ? refMultHalfPosWest++ : refMultHalfNegWest++;
            }
            if (pMom.Eta() < 0) {
               (charge > 0) ? refMultHalfPosEast++ : refMultHalfNegEast++;
            }

            if (0.5 < pMom.Eta() && pMom.Eta() < 1.) {
               (charge > 0) ? refMult2PosWest++ : refMult2NegWest++;
            }
            if (-1. < pMom.Eta() && pMom.Eta() < -0.5) {
               (charge > 0) ? refMult2PosEast++ : refMult2NegEast++;
            }

         } // if ( miniTrack->nHits() > 15 )
      }    // if (isPrimary)
      else {
         miniTrack->setPrimaryMomentum(TVector3(0., 0., 0.));
      }

      // Retrieve global track to get more parameters that are not available
      // directly from kalman track
      MpdTrack *glTrack = (MpdTrack *)glTracks->UncheckedAt(iKalmanTrk);

      // Track must exist
      if (!glTrack || (kalTrack->GetTrackID() != glTrack->GetID())) {
         std::cout << "[WARNING] MpdMiniDstFillTask::fillTracks -"
                   << " No gtrk that corresponds to ktrk" << std::endl;
         continue;
      }

      // Setting global track momentum at DCA to primary vertex
      TVector3 globMom(glTrack->GetPx(), glTrack->GetPy(), glTrack->GetPz());
      miniTrack->setGlobalMomentum(globMom.X(), globMom.Y(), globMom.Z());

      // Set nSigma and dE/dx info
      if (fNSigmaDedxEstimator == 0) {
         miniTrack->setNSigmaElectron(glTrack->GetNSigmaElectron());
         miniTrack->setNSigmaPion(glTrack->GetNSigmaPion());
         miniTrack->setNSigmaKaon(glTrack->GetNSigmaKaon());
         miniTrack->setNSigmaProton(glTrack->GetNSigmaProton());
      } else if (fNSigmaDedxEstimator == 1) {
         nSigma = nSigmaDedx(globMom.Mag(), glTrack->GetdEdXTPC());
         miniTrack->setNSigmaElectron(nSigma.at(0));
         miniTrack->setNSigmaPion(nSigma.at(1));
         miniTrack->setNSigmaKaon(nSigma.at(2));
         miniTrack->setNSigmaProton(nSigma.at(3));
      } else {
         std::cout << "[WARNING] Wrong fNSigmaDedxEstimator: " << fNSigmaDedxEstimator << std::endl;
      }
      miniTrack->setDedx(glTrack->GetdEdXTPC());
      miniTrack->setHitMap(glTrack->GetLayerHitMap());

      // Getting primary vertex the first primary vertex
      MpdVertex *vertex = (MpdVertex *)fVertices->UncheckedAt(0);

      // Get primary vertex position
      TVector3 primVertex(vertex->GetX(), vertex->GetY(), vertex->GetZ());
      TVector3 firstPoint(glTrack->GetFirstPointX(), glTrack->GetFirstPointY(), glTrack->GetFirstPointZ());

      // Physical helix instantiation
      MpdMiniPhysicalHelix helix(globMom, firstPoint, fBField * kilogauss, glTrack->GetCharge());
      double               pathLength  = helix.pathLength(primVertex);
      TVector3             dcaPosition = helix.at(pathLength);
      miniTrack->setOrigin(dcaPosition.X(), dcaPosition.Y(), dcaPosition.Z());

      // Estimate refMult of global tracks
      if (TMath::Abs(miniTrack->gMom().Eta()) < 0.5 && miniTrack->charge() != 0 && miniTrack->nHits() > 15) {
         grefMult++;
      }

      // Store index to miniMcTrack
      miniTrack->setMcTrackIndex(miniMcIdxFromMcIdx(kalTrack->GetTrackID()));

      // Update miniMcTrack global track ID information
      // in case of match
      if (miniTrack->mcTrackIndex() >= 0) {

         MpdMiniMcTrack *miniMcTrack =
            (MpdMiniMcTrack *)fMiniArrays[MpdMiniArrays::McTrack]->At(miniMcIdxFromMcIdx(kalTrack->GetTrackID()));
         miniMcTrack->addGlobalTrackId((UShort_t)miniTracks->GetEntriesFast() - 1);
      }

      //
      //
      // ECal-matching information
      //
      //

      auto itr2 = fMcTrk2EcalCluster.find(kalTrack->GetTrackID());
      if (itr2 != fMcTrk2EcalCluster.end()) {
         miniTrack->setBECalClusterIndex(itr2->second);
      } else {
         miniTrack->setBECalClusterIndex(-1);
      }

      // Check number of BECal-matched tracks
      if (miniTrack->isBECalTrack()) nBECalMatched++;

      //
      //
      // TOF-matching information
      //
      // TODO: must recalculate the information for primary tracks
      //       because it does not make any sense for kalman or global ones
      //
      //

      // Loop over TOF-matching information
      if (fTofMatching)
         for (Int_t iMatch = 0; iMatch < fTofMatching->GetEntriesFast(); iMatch++) {

            // Retrieve TOF-matching information
            MpdTofMatchingData *dataMatch = (MpdTofMatchingData *)fTofMatching->UncheckedAt(iMatch);

            // TOF matching information should exist
            if (!dataMatch || (dataMatch->GetKFTrackIndex() != iKalmanTrk)) continue;

            // Create barrel TOF traits
            MpdMiniBTofPidTraits *bTofPidTraits =
               new ((*miniBTofTraits)[miniBTofTraits->GetEntriesFast()]) MpdMiniBTofPidTraits();

            miniTrack->setBTofPidTraitsIndex(miniBTofTraits->GetEntriesFast() - 1);

            bTofPidTraits->setTrackIndex(miniTracks->GetEntriesFast() - 1);
            bTofPidTraits->setHitIndex(dataMatch->GetTofHitIndex());
            bTofPidTraits->setBeta(dataMatch->GetBeta());
            bTofPidTraits->setMomentum(dataMatch->GetMomentum());

            // Increment number of TOF-matched tracks
            nTofMatched++;

            break;
         } // for (Int_t iMatch = 0; iMatch < fTofMatching->GetEntriesFast(); iMatch++)

   } // for (Int_t iKalmanTrk = 0; iKalmanTrk < fTpcTracks->GetEntriesFast(); iKalmanTrk++)

   //
   // Fill MpdMiniEvent with reference multiplicity values
   //
   MpdMiniEvent *miniEvent = (MpdMiniEvent *)fMiniArrays[MpdMiniArrays::Event]->At(0);
   miniEvent->setRefMultNeg(refMultNeg);
   miniEvent->setRefMultPos(refMultPos);
   miniEvent->setRefMultHalfNegEast(refMultHalfNegEast);
   miniEvent->setRefMultHalfPosEast(refMultHalfPosEast);
   miniEvent->setRefMultHalfNegWest(refMultHalfNegWest);
   miniEvent->setRefMultHalfPosWest(refMultHalfPosWest);
   miniEvent->setRefMult2NegEast(refMult2NegEast);
   miniEvent->setRefMult2PosEast(refMult2PosEast);
   miniEvent->setRefMult2NegWest(refMult2NegWest);
   miniEvent->setRefMult2PosWest(refMult2PosWest);
   miniEvent->setGRefMult(grefMult);
   miniEvent->setNumberOfBTOFMatch(nTofMatched);
   miniEvent->setNumberOfBECalMatch(nBECalMatched);
}

//_________________

void MpdMiniDstFillTask::refit2Vp(MpdMiniTrack *miniTrack, Int_t iTpcKalmanTrack, MpdVertex *vtx)
{

   // Get primary tracks from event
   // Done by smooth tracks from primary vertex (update momentum and track length -
   // covariance matrix is not updated !!!)
   // Got from MpdKfPrimaryVertexFinder::Smooth()
   MpdKalmanHit hit;
   TMatrixD     c(3, 3), xk(3, 1), ck0(5, 1);
   TMatrixD     a(5, 3), b(5, 3);
   TVector3     vert;
   vtx->Position(vert);
   xk(0, 0)     = vert.X();
   xk(1, 0)     = vert.Y();
   xk(2, 0)     = vert.Z();
   Double_t rad = vert.Pt();

   MpdKalmanTrack *track  = (MpdKalmanTrack *)fTpcTracks->UncheckedAt(iTpcKalmanTrack);
   MpdKalmanTrack *trVert = new MpdKalmanTrack(*track);

   MpdKalmanTrack track1 = *track;
   track1.SetParamNew(*track1.GetParam());
   track1.SetPos(track1.GetPosNew());
   track1.ReSetWeight();
   track1.SetLength(0.);
   TMatrixD g = *track1.GetWeight(); // track weight matrix

   if (track->GetNode() == "") {
      hit.SetType(MpdKalmanHit::kFixedR);
      hit.SetPos(track->GetPos());
   } else {
      hit.SetType(MpdKalmanHit::kFixedP);
      TString detName = track->GetNode();
      if (track->GetUniqueID()) {
         // ITS
         detName = detName(16, detName.Length());
         detName += "#0";
      }
      MpdKalmanGeoScheme *geo = MpdKalmanFilter::Instance()->GetGeo();
      hit.SetDetectorID(geo->DetId(detName));

      // Find distance from the current track position to the last point (plane) -
      // to define direction (mainly for ITS)
      TVector3 pos   = geo->GlobalPos(&hit);
      TVector3 norm  = geo->Normal(&hit);
      Double_t v7[7] = {0.0};
      track1.SetNode("");
      MpdKalmanFilter::Instance()->SetGeantParamB(&track1, v7, 1);
      Double_t d = -(pos * norm); // Ax+By+Cz+D=0, A=nx, B=ny, C=nz
      TVector3 v3(v7[0], v7[1], v7[2]);
      d += v3 * norm;
      if (d < 0) {
         track1.SetDirection(MpdKalmanTrack::kOutward);
      }
   } // else

   MpdKalmanFilter::Instance()->PropagateToHit(&track1, &hit, kTRUE, kTRUE);

   // TMatrixD* par2 = track->GetParamNew();
   // cout << par2 << endl;
   computeAandB(xk, track, track1, a, b, ck0); // compute matrices of derivatives

   // W = (Bt*G*B)'
   TMatrixD tmp(g, TMatrixD::kMult, b);
   TMatrixD w(b, TMatrixD::kTransposeMult, tmp);
   w.Invert();

   TMatrixD m = *track1.GetParamNew();
   m -= ck0; // m-ck0

   // qk = W*Bt*G*(m-ck0-A*xk)
   TMatrixD tmp21(a, TMatrixD::kMult, xk);
   tmp21 *= -1;
   tmp21 += m; // m-ck0-A*xk
   TMatrixD tmp22(g, TMatrixD::kMult, tmp21);
   TMatrixD tmp23(b, TMatrixD::kTransposeMult, tmp22);
   TMatrixD qk(w, TMatrixD::kMult, tmp23);

   // Update momentum and last coordinate
   TMatrixD *parPointer = nullptr;

   if (track->GetParamNew())
      parPointer = track->GetParamNew();
   else
      parPointer = &m;

   TMatrixD par = *parPointer;
   for (Int_t i = 0; i < 3; ++i) par(i + 2, 0) = qk(i, 0);
   par(0, 0) = rad * vert.Phi();
   par(1, 0) = vert.Z();
   trVert->SetParam(par);
   trVert->SetPosNew(rad);

   Double_t pT    = 1. / trVert->GetParam(4); // Signed pT
   Double_t phi   = trVert->GetParam(2);
   Double_t theta = TMath::PiOver2() - trVert->GetParam(3);

   Double_t Px = TMath::Abs(pT) * TMath::Cos(phi);
   Double_t Py = TMath::Abs(pT) * TMath::Sin(phi);

   Double_t Pz = -1.;
   if (TMath::Sin(theta) != 0.) {
      Pz = TMath::Abs(pT) / TMath::Tan(theta);
   }

   miniTrack->setPrimaryMomentum(TVector3(Px, Py, Pz));
}

//_________________

void MpdMiniDstFillTask::computeAandB(TMatrixD &xk0, const MpdKalmanTrack *track, const MpdKalmanTrack &trackM,
                                      TMatrixD &a, TMatrixD &b, TMatrixD &ck0)
{

   // Compute matrices of derivatives w.r.t. vertex coordinates and track momentum
   Double_t vert0[3], *vert = xk0.GetMatrixArray(); // zero[3] = {0},
   for (Int_t i = 0; i < 3; ++i) vert0[i] = vert[i];

   MpdKalmanTrack trackk = *track;
   trackk.SetPos(trackk.GetPosNew());
   // trackk.GetParam()->Print();
   // Propagate track to PCA w.r.t. point xk0
   MpdKalmanFilter::Instance()->FindPca(&trackk, vert0);
   // MpdKalmanFilter::Instance()->FindPca(&trackk,zero); // just for test
   // std::cout << trackk.GetPosNew() << std::endl;
   trackk.SetParam(*trackk.GetParamNew());
   // trackk.GetParam()->Print();

   // Put track at xk0
   Double_t r   = TMath::Sqrt(vert0[0] * vert0[0] + vert0[1] * vert0[1]);
   Double_t phi = trackk.GetParamNew(2); // track Phi
   if (r > 1.e-7) phi = TMath::ATan2(vert0[1], vert0[0]);
   trackk.SetPos(r);
   trackk.SetParam(0, r * phi);
   trackk.SetParam(1, vert0[2]);
   trackk.SetNode("");
   MpdKalmanTrack track0 = trackk;

   // Propagate track to chosen radius
   MpdKalmanHit hit;
   // hit.SetR(35.);
   // hit = *(MpdKalmanHitR*)track->GetTrHits()->Last();
   if (track->GetNode() == "") {
      hit.SetType(MpdKalmanHit::kFixedR);
      // hit.SetR(35.);
      // hit = *(MpdKalmanHitR*)track->GetTrHits()->Last();
      hit.SetPos(track->GetPos());
      MpdKalmanFilter::Instance()->PropagateParamR(&trackk, &hit, kFALSE);
      // trackk.GetParamNew()->Print();
      Proxim(trackM, trackk);
   } else {
      hit.SetType(MpdKalmanHit::kFixedP);
      TString detName = track->GetNode();
      if (track->GetUniqueID()) {
         // ITS
         detName = detName(16, detName.Length());
         detName += "#0";
      }
      MpdKalmanGeoScheme *geo = MpdKalmanFilter::Instance()->GetGeo();
      hit.SetDetectorID(geo->DetId(detName));
      // Find distance from the current track position to the last point (plane) -
      // to define direction (mainly for ITS)
      TVector3 pos   = geo->GlobalPos(&hit);
      TVector3 norm  = geo->Normal(&hit);
      Double_t v7[7] = {0.0};
      MpdKalmanFilter::Instance()->SetGeantParamB(&trackk, v7, 1);
      Double_t d = -(pos * norm); // Ax+By+Cz+D=0, A=nx, B=ny, C=nz
      TVector3 v3(v7[0], v7[1], v7[2]);
      d += v3 * norm;
      if (d < 0) trackk.SetDirection(MpdKalmanTrack::kOutward);
      MpdKalmanFilter::Instance()->PropagateParamP(&trackk, &hit, kFALSE, kTRUE);
      track0.SetDirection(trackk.GetDirection());
   }

   // Double_t shift = 0.01; // 100 um coordinate shift
   Double_t shift = 0.1; // 1 mm coordinate shift
   for (Int_t i = 0; i < 3; ++i) {
      MpdKalmanTrack track1 = track0;
      vert0[i] += shift;
      if (i > 0) vert0[i - 1] -= shift;
      r = TMath::Sqrt(vert0[0] * vert0[0] + vert0[1] * vert0[1]);
      if (r > 1.e-7)
         phi = TMath::ATan2(vert0[1], vert0[0]);
      else
         phi = track0.GetParamNew(2); // track Phi
      track1.SetPos(r);
      track1.SetParam(0, r * phi);
      track1.SetParam(1, vert0[2]);
      if (track->GetNode() == "") {
         MpdKalmanFilter::Instance()->PropagateParamR(&track1, &hit, kFALSE);
         Proxim(trackk, track1);
         // proxim(track1,trackk);
      } else
         MpdKalmanFilter::Instance()->PropagateParamP(&track1, &hit, kFALSE, kTRUE);
      // Derivatives
      for (Int_t j = 0; j < 5; ++j) {
         a(j, i) = (track1.GetParamNew(j) - trackk.GetParamNew(j)) / shift;
      }
   } // for (Int_t i = 0; i < 3; ++i)

   for (Int_t i = 0; i < 3; ++i) {
      MpdKalmanTrack track1 = track0;
      Int_t          j      = i + 2;
      shift                 = (*track->GetCovariance())(j, j);
      shift                 = TMath::Sqrt(shift);
      if (j == 4) shift *= TMath::Sign(1., -track0.GetParamNew(j)); // 1/p
      track1.SetParam(j, track0.GetParamNew(j) + shift);
      // if (j == 2 && track1.GetParamNew(j)*TMath::Sign(1.,track1.GetParamNew(j)) > TMath::Pi())
      // track1.SetParam(j,track0.GetParamNew(j)-shift);
      if (track->GetNode() == "") {
         MpdKalmanFilter::Instance()->PropagateParamR(&track1, &hit, kFALSE);
         Proxim(trackk, track1);
         // proxim(track1,trackk);
      } else
         MpdKalmanFilter::Instance()->PropagateParamP(&track1, &hit, kFALSE, kTRUE);
      // Derivatives
      for (Int_t k = 0; k < 5; ++k) {
         b(k, i) = (track1.GetParamNew(k) - trackk.GetParamNew(k)) / shift;
      }
   } // for (Int_t i = 0; i < 3; ++i)

   TMatrixD qk0(3, 1);
   for (Int_t i = 0; i < 3; ++i) qk0(i, 0) = track0.GetParamNew(i + 2);
   // qk0.Print();
   ck0 = *trackk.GetParamNew();
   ck0 -= TMatrixD(a, TMatrixD::kMult, xk0);
   ck0 -= TMatrixD(b, TMatrixD::kMult, qk0);
}

//________________

void MpdMiniDstFillTask::Proxim(const MpdKalmanTrack &track0, MpdKalmanTrack &track)
{

   if (track0.GetType() != MpdKalmanTrack::kBarrel) {
      std::cout << " !!! Implemented only for kBarrel tracks !!!" << std::endl;
      exit(0);
   }

   // Double_t tmp = track.GetParamNew(0);
   Double_t phi0 = track0.GetParamNew(0) / track0.GetPosNew();
   Double_t phi  = track.GetParamNew(0) / track.GetPosNew();
   phi           = MpdKalmanFilter::Instance()->Proxim(phi0, phi);
   TMatrixD *par = track.GetParamNew();
   (*par)(0, 0)  = phi * track.GetPosNew();
   phi0          = track0.GetParamNew(2);
   phi           = track.GetParamNew(2);
   phi           = MpdKalmanFilter::Instance()->Proxim(phi0, phi);
   (*par)(2, 0)  = phi;
   track.SetParamNew(*par);
}

//_________________

void MpdMiniDstFillTask::fillCovMatrix(MpdTpcKalmanTrack *tpcTrack, MpdMiniTrackCovMatrix *miniTrackCovMatrix)
{

   // Covariance matrix is stored at DCA2beamline as being written in the source track !!!

   Double_t rPhi0    = tpcTrack->GetParam()->GetMatrixArray()[0];
   Double_t Z        = tpcTrack->GetParam()->GetMatrixArray()[1];
   Double_t Phi      = tpcTrack->GetParam()->GetMatrixArray()[2];
   Double_t Lambda   = tpcTrack->GetParam()->GetMatrixArray()[3];
   Double_t MinusQPt = tpcTrack->GetParam()->GetMatrixArray()[4];

   Double_t R  = tpcTrack->GetPosNew();
   Double_t R0 = tpcTrack->GetPos();

   miniTrackCovMatrix->SetStateVector(rPhi0, Z, Phi, Lambda, MinusQPt);
   miniTrackCovMatrix->SetR(R);
   miniTrackCovMatrix->SetR0(R0);

   vector<Double_t> sigmas;
   vector<Double_t> correlations;

   for (Int_t i = 0; i < 5; i++)
      for (Int_t j = 0; j < 5; j++) {
         if (i == j)
            sigmas.push_back((*tpcTrack->GetCovariance())[i][j]);
         else if (i > j)
            correlations.push_back((*tpcTrack->GetCovariance())[i][j]);
      }

   miniTrackCovMatrix->setSigmas(sigmas);
   miniTrackCovMatrix->setCorrelations(correlations);
}

//_________________

Int_t MpdMiniDstFillTask::miniMcIdxFromMcIdx(Int_t mcIdx)
{
   Int_t ret = -1; // Not found
   for (UInt_t i = 0; i < fMcTrk2MiniMcTrk.size(); i++) {
      if (fMcTrk2MiniMcTrk.at(i).first == mcIdx) {
         ret = fMcTrk2MiniMcTrk.at(i).second;
         break;
      }
   }
   return ret;
}

//________________

std::vector<Double_t> MpdMiniDstFillTask::nSigmaDedx(Double_t p, Double_t dEdx)
{

   std::vector<Float_t> electronMean{
      3757.34, 3770.90, 3784.94, 3799.03, 3808.90, 3819.55, 3828.35, 3839.03, 3851.42, 3860.54, 3868.62, 3875.06,
      3883.33, 3891.81, 3898.33, 3906.15, 3912.76, 3920.30, 3924.64, 3929.60, 3935.65, 3939.86, 3942.09, 3949.23,
      3953.97, 3958.03, 3961.38, 3966.54, 3967.48, 3972.15, 3974.84, 3976.74, 3982.19, 3983.85, 3989.49, 3990.05,
      3990.49, 3996.64, 3996.57, 3998.25, 4000.55, 4004.68, 4004.29, 4009.05, 4011.62, 4018.04, 4016.26, 4016.22,
      4019.27, 4017.20, 4020.18, 4024.86, 4023.67, 4025.07, 4028.15, 4031.01, 4032.08, 4035.21, 4029.78, 4043.04,
      4039.36, 4040.60, 4042.16, 4041.30, 4034.07, 4046.29, 4043.97, 4052.67, 4050.12, 4055.28, 4053.21, 4047.72,
      4051.43, 4050.56, 4057.61, 4061.93, 4053.38, 4062.25, 4066.25, 4063.20, 4064.63, 4066.04, 4066.93, 4074.67,
      4080.85, 4067.16, 4059.44, 4065.39, 4075.59, 4078.03, 4071.79, 4072.62, 4085.47, 4077.24, 4085.17, 4085.90,
      4067.84, 4085.41, 4094.88, 4084.99, 4079.28, 4090.99, 4085.75, 4109.16, 4085.81};
   std::vector<Float_t> pionMean{
      4262.84, 4036.03, 3845.60, 3719.01, 3577.86, 3465.71, 3393.91, 3310.76, 3226.63, 3177.69, 3136.60, 3089.69,
      3041.34, 3000.81, 2984.68, 2988.91, 3006.90, 2988.47, 2963.48, 2946.04, 2939.60, 2928.94, 2917.17, 2905.05,
      2892.70, 2880.60, 2870.70, 2868.63, 2866.40, 2862.11, 2857.59, 2852.92, 2848.26, 2843.56, 2839.01, 2835.19,
      2834.63, 2836.07, 2835.99, 2835.48, 2834.87, 2834.14, 2833.61, 2833.06, 2832.32, 2831.84, 2832.19, 2833.61,
      2834.33, 2835.11, 2835.17, 2835.92, 2836.62, 2837.92, 2839.48, 2840.95, 2842.79, 2844.75, 2846.28, 2848.05,
      2849.51, 2851.45, 2853.13, 2854.91, 2856.73, 2858.61, 2859.91, 2861.68, 2863.49, 2865.84, 2867.34, 2869.98,
      2871.76, 2873.91, 2876.03, 2877.96, 2880.20, 2882.43, 2884.52, 2886.60, 2888.32, 2890.91, 2892.97, 2895.19,
      2897.32, 2899.75, 2901.64, 2903.70, 2906.02, 2908.15, 2910.38, 2912.85, 2914.53, 2916.71, 2918.79, 2920.87,
      2922.91, 2925.07, 2928.05, 2929.70, 2931.79, 2933.62, 2936.53, 2938.11, 2940.14, 2942.70, 2944.81, 2947.09,
      2948.71, 2951.36, 2953.11, 2955.20, 2957.73, 2959.67, 2962.50, 2963.49, 2965.71, 2967.79, 2969.70, 2971.68,
      2974.71, 2976.76, 2978.99, 2980.92, 2982.49, 2983.69, 2986.54, 2987.96, 2989.94, 2992.44, 2994.25, 2996.08,
      2996.84, 3000.75, 3002.46, 3002.93, 3005.36, 3008.21, 3009.69, 3011.52, 3012.96, 3014.92, 3016.88, 3018.95,
      3021.46, 3023.63, 3025.01, 3027.19, 3027.37, 3030.46, 3032.71, 3033.61, 3037.03, 3037.25, 3039.41, 3041.14,
      3042.50, 3046.58, 3049.13, 3049.04, 3051.06, 3052.68, 3053.95, 3054.28, 3057.88, 3057.84, 3059.71, 3060.80,
      3063.40, 3067.77, 3066.87, 3069.38, 3071.58, 3071.18, 3074.22, 3074.04, 3077.39, 3077.13, 3079.85, 3081.78,
      3082.59, 3085.07, 3086.51, 3088.73, 3089.27, 3089.88, 3093.16, 3095.17, 3094.81, 3097.77, 3096.85, 3102.71,
      3100.34, 3103.17, 3105.11, 3108.30, 3107.99, 3110.87, 3112.43, 3112.21, 3113.52, 3116.44, 3118.44, 3122.84,
      3119.92, 3124.31, 3125.29, 3124.75, 3125.22, 3128.23, 3128.21, 3126.40, 3128.76, 3131.73, 3131.05, 3135.59,
      3136.53, 3133.55, 3137.90, 3136.49, 3140.83, 3143.74, 3138.88, 3142.34, 3147.02, 3140.30, 3138.59, 3139.91,
      3139.16, 3143.03, 3146.28, 3137.80, 3134.58, 3137.14, 3138.64};
   std::vector<Float_t> kaonMean{
      24902.63, 22512.42, 20653.74, 18948.00, 17529.94, 16309.94, 15078.43, 14133.55, 13229.39, 12313.20, 11611.54,
      11000.31, 10330.84, 9745.07,  9322.95,  8911.44,  8465.27,  8038.40,  7739.63,  7482.00,  7201.13,  6908.41,
      6630.43,  6425.87,  6264.29,  6089.32,  5904.01,  5713.26,  5533.83,  5404.25,  5301.87,  5191.32,  5074.85,
      4955.59,  4838.39,  4726.18,  4637.42,  4578.42,  4512.20,  4442.14,  4369.90,  4295.83,  4220.43,  4145.17,
      4078.87,  4030.65,  3993.99,  3955.50,  3913.45,  3870.88,  3825.12,  3779.50,  3734.01,  3689.94,  3645.96,
      3611.15,  3585.50,  3565.63,  3543.24,  3518.78,  3493.96,  3468.58,  3441.34,  3415.68,  3388.67,  3362.41,
      3336.47,  3312.95,  3293.65,  3278.79,  3267.47,  3253.99,  3240.52,  3227.13,  3214.76,  3198.35,  3184.53,
      3170.68,  3156.23,  3142.34,  3128.98,  3116.82,  3106.47,  3096.47,  3085.51,  3074.93,  3066.53,  3059.45,
      3054.90,  3053.29,  3059.15,  3063.03,  3066.18,  3068.04,  3066.32,  3062.43,  3057.32,  3054.12,  3047.56,
      3042.65,  3038.70,  3034.49,  3028.39,  3023.09,  3016.85,  3013.23,  3006.28,  3000.60,  2996.01,  2990.71,
      2986.44,  2981.70,  2978.35,  2973.94,  2969.85,  2966.46,  2963.08,  2963.20,  2960.06,  2956.58,  2952.46,
      2949.20,  2949.07,  2944.48,  2943.52,  2941.92,  2939.51,  2935.96,  2934.03,  2933.08,  2929.91,  2927.59,
      2925.62,  2925.12,  2920.79,  2918.91,  2917.29,  2916.21,  2914.11,  2909.00,  2910.72,  2906.60,  2907.18,
      2906.63,  2903.04,  2902.20,  2903.17,  2902.20,  2902.07,  2898.14,  2896.82,  2899.09,  2897.26,  2898.67,
      2895.51,  2896.22,  2895.20,  2894.95,  2896.60,  2892.48,  2894.59,  2893.13,  2891.10,  2890.79,  2892.82,
      2889.15,  2888.74,  2888.66,  2885.85,  2883.53,  2886.15,  2885.93,  2885.98,  2885.17,  2883.49,  2884.08,
      2883.31,  2884.29,  2880.32,  2885.39,  2885.55,  2884.25,  2881.87,  2883.70,  2883.74,  2884.78,  2884.68,
      2887.26,  2888.48,  2881.29,  2886.50,  2882.68,  2880.41,  2887.89,  2885.83,  2887.97,  2885.93,  2887.57,
      2884.88,  2883.68,  2889.82,  2882.39,  2884.60,  2883.05,  2880.41,  2887.52,  2885.94,  2883.64,  2879.01,
      2886.20,  2887.72,  2886.49,  2888.96,  2880.14,  2883.86,  2889.81,  2890.89,  2888.92,  2900.42,  2890.46,
      2892.30,  2893.44,  2885.62,  2895.89,  2890.43,  2891.32,  2895.11,  2890.97,  2889.99,  2898.71,  2892.84,
      2896.20,  2889.65,  2897.49,  2887.01};
   std::vector<Float_t> protonMean{
      38972.37, 37573.08, 35969.80, 34459.91, 33019.80, 31650.71, 30150.28, 28471.70, 26478.74, 24943.29, 23634.93,
      22394.18, 21419.11, 20571.20, 19735.88, 18887.98, 18096.00, 17437.12, 16846.17, 16207.97, 15554.46, 14939.65,
      14434.02, 14004.97, 13563.70, 13095.23, 12612.15, 12169.34, 11803.16, 11504.49, 11193.70, 10861.69, 10516.00,
      10173.81, 9866.30,  9626.29,  9428.65,  9223.85,  9006.58,  8783.11,  8547.25,  8314.05,  8099.92,  7924.99,
      7789.16,  7658.45,  7516.56,  7369.65,  7219.12,  7064.38,  6911.46,  6761.07,  6628.94,  6523.83,  6437.91,
      6351.78,  6261.61,  6166.89,  6070.43,  5970.76,  5870.26,  5769.79,  5671.57,  5581.27,  5506.87,  5446.44,
      5392.92,  5336.64,  5277.38,  5217.60,  5155.84,  5093.41,  5030.03,  4966.69,  4905.62,  4842.56,  4783.99,
      4734.08,  4691.10,  4654.32,  4622.63,  4587.78,  4552.15,  4516.16,  4478.47,  4439.35,  4400.60,  4361.30,
      4321.57,  4282.49,  4243.61,  4204.60,  4167.69,  4133.43,  4104.35,  4078.59,  4057.82,  4037.31,  4017.78,
      3996.53,  3973.57,  3952.16,  3929.09,  3905.43,  3882.12,  3858.92,  3836.02,  3812.96,  3790.78,  3769.17,
      3749.01,  3728.11,  3708.20,  3689.55,  3670.49,  3652.33,  3634.04,  3615.72,  3597.41,  3581.73,  3564.05,
      3548.06,  3533.77,  3519.34,  3505.96,  3494.08,  3481.90,  3471.18,  3460.27,  3448.81,  3437.92,  3427.87,
      3417.28,  3406.27,  3395.20,  3386.34,  3375.71,  3365.49,  3355.56,  3344.80,  3334.77,  3323.59,  3314.08,
      3303.63,  3294.71,  3285.58,  3275.74,  3265.11,  3256.75,  3249.46,  3241.93,  3235.16,  3227.95,  3219.91,
      3215.56,  3208.66,  3202.83,  3197.98,  3190.52,  3185.80,  3180.11,  3174.12,  3168.97,  3162.83,  3157.70,
      3152.17,  3146.97,  3141.88,  3137.22,  3131.89,  3126.64,  3122.98,  3117.68,  3114.53,  3109.62,  3106.55,
      3104.21,  3102.82,  3101.72,  3102.88,  3101.63,  3103.90,  3104.95,  3105.39,  3106.98,  3107.46,  3107.23,
      3107.80,  3107.03,  3105.42,  3105.41,  3105.67,  3101.84,  3099.91,  3098.73,  3094.83,  3093.52,  3091.87,
      3088.97,  3086.46,  3082.45,  3080.88,  3079.28,  3075.34,  3072.60,  3071.95,  3068.22,  3066.21,  3063.35,
      3061.61,  3057.89,  3056.60,  3052.09,  3051.42,  3048.96,  3045.14,  3040.91,  3040.11,  3035.33,  3035.78,
      3032.39,  3030.79,  3028.90,  3028.71,  3024.20,  3024.21,  3020.38,  3018.97,  3015.58,  3015.54,  3010.61};
   std::vector<Float_t> electronSigma{
      258.53, 258.66, 258.60, 259.82, 260.60, 262.01, 262.63, 264.18, 266.58, 265.97, 266.51, 266.67, 267.31, 268.45,
      266.98, 269.08, 269.89, 270.52, 268.64, 269.40, 271.10, 270.64, 271.42, 272.45, 270.73, 273.04, 274.42, 275.20,
      276.39, 276.95, 277.20, 278.11, 278.62, 279.04, 279.23, 279.30, 279.41, 278.45, 280.47, 282.71, 282.09, 282.52,
      282.94, 283.31, 286.75, 285.68, 289.90, 285.61, 286.83, 284.81, 288.44, 285.30, 286.75, 288.75, 289.52, 287.85,
      288.96, 293.38, 292.45, 293.87, 291.02, 293.36, 294.46, 290.53, 291.99, 297.53, 297.21, 293.31, 291.58, 296.65,
      301.40, 294.78, 294.42, 299.74, 295.31, 297.58, 293.69, 300.13, 293.72, 298.41, 289.15, 299.23, 311.42, 297.80,
      305.38, 304.40, 301.76, 296.18, 302.25, 290.06, 307.85, 302.16, 313.76, 303.03, 298.38, 307.94, 297.13, 303.07,
      296.57, 290.59, 303.75, 293.38, 300.15, 292.96, 289.06};
   std::vector<Float_t> pionSigma{
      343.49, 323.76, 291.77, 277.81, 262.79, 244.54, 237.33, 229.45, 220.80, 213.92, 210.48, 206.63, 202.53, 197.83,
      195.85, 194.98, 194.86, 193.02, 191.03, 189.22, 188.77, 188.25, 187.49, 186.59, 185.96, 185.08, 184.25, 184.19,
      184.34, 183.99, 183.92, 183.57, 183.26, 183.17, 182.81, 182.70, 182.59, 182.87, 183.16, 183.00, 183.06, 183.28,
      183.26, 183.44, 183.39, 183.57, 183.65, 183.85, 184.11, 184.31, 184.15, 184.38, 184.55, 184.78, 185.25, 185.43,
      185.87, 186.04, 186.06, 186.19, 186.61, 186.95, 187.17, 187.46, 187.63, 187.92, 187.99, 188.42, 188.61, 188.81,
      189.10, 189.43, 189.81, 189.93, 190.42, 190.40, 190.58, 191.20, 191.00, 191.89, 191.55, 191.65, 192.39, 192.60,
      192.70, 192.92, 193.21, 194.02, 193.96, 193.80, 194.05, 194.42, 195.15, 194.74, 195.15, 195.86, 195.69, 196.13,
      196.48, 196.64, 197.03, 197.45, 197.77, 197.26, 197.79, 197.93, 197.78, 197.97, 198.55, 198.90, 199.04, 199.16,
      199.54, 199.16, 199.66, 200.10, 200.39, 200.81, 200.45, 200.24, 200.93, 201.21, 201.48, 201.68, 201.68, 201.61,
      202.66, 202.52, 202.05, 201.35, 201.89, 202.37, 203.61, 203.45, 203.75, 203.88, 203.58, 203.09, 203.31, 204.57,
      204.46, 204.89, 205.14, 205.81, 205.79, 204.50, 206.12, 206.07, 205.86, 206.40, 206.25, 206.56, 206.88, 206.72,
      206.06, 206.89, 208.78, 206.17, 208.23, 206.99, 208.42, 208.33, 207.79, 207.94, 207.77, 206.79, 208.69, 209.45,
      209.40, 211.47, 209.18, 209.05, 211.09, 208.65, 210.92, 211.11, 210.44, 211.16, 210.46, 211.39, 210.92, 210.99,
      210.44, 211.89, 212.97, 211.83, 210.25, 211.96, 212.69, 213.38, 212.21, 215.57, 213.26, 214.98, 215.59, 214.87,
      216.26, 216.37, 215.11, 216.17, 216.37, 213.79, 216.14, 214.42, 215.22, 214.40, 215.69, 215.58, 217.22, 216.90,
      219.18, 217.05, 217.72, 220.35, 220.81, 219.64, 221.73, 220.45, 219.60, 220.21, 220.55, 216.99, 226.87, 220.02,
      223.34, 218.94, 221.12, 225.28, 226.55, 224.69, 226.13, 224.44, 221.95, 226.09, 226.29};
   std::vector<Float_t> kaonSigma{
      3134.25, 2584.22, 2206.67, 1875.84, 1560.24, 1461.16, 1237.64, 1089.92, 1079.11, 966.25, 825.77, 823.58, 805.53,
      699.73,  638.04,  639.51,  631.65,  575.04,  518.67,  510.84,  506.91,  497.52,  463.83, 427.75, 419.60, 417.30,
      412.28,  402.13,  383.97,  359.91,  353.49,  350.35,  347.75,  342.51,  334.93,  324.09, 309.56, 304.00, 304.51,
      302.02,  299.22,  296.38,  293.62,  286.25,  278.27,  271.87,  268.80,  268.85,  266.46, 264.79, 263.43, 261.02,
      259.30,  254.83,  251.58,  246.21,  243.71,  242.35,  241.68,  240.91,  240.07,  238.66, 237.05, 236.49, 234.23,
      232.20,  229.64,  227.52,  224.01,  223.42,  222.02,  222.10,  222.17,  220.40,  220.11, 219.90, 219.86, 217.23,
      216.55,  215.35,  214.89,  213.22,  212.78,  210.81,  209.93,  209.11,  208.96,  207.31, 206.47, 205.72, 207.58,
      207.22,  207.81,  208.01,  206.87,  206.50,  206.42,  206.03,  205.53,  205.32,  204.97, 204.17, 203.76, 204.34,
      202.50,  203.19,  202.43,  201.42,  201.10,  201.23,  201.14,  201.00,  199.74,  200.99, 197.84, 199.02, 198.92,
      198.01,  198.98,  198.91,  198.03,  198.26,  199.48,  198.51,  198.77,  197.56,  199.30, 197.04, 196.83, 197.60,
      197.61,  196.84,  196.54,  197.20,  196.71,  196.54,  195.94,  196.31,  194.34,  195.07, 196.95, 195.58, 194.33,
      195.75,  194.91,  196.49,  196.18,  196.57,  195.60,  195.58,  195.65,  197.12,  196.27, 195.39, 196.52, 194.10,
      195.39,  195.90,  194.37,  197.73,  195.24,  198.60,  194.33,  193.05,  193.43,  195.81, 193.50, 196.00, 197.59,
      194.61,  196.57,  196.15,  197.28,  192.05,  192.17,  195.60,  195.71,  194.63,  198.89, 195.48, 193.99, 194.54,
      196.01,  197.04,  198.56,  193.06,  200.18,  193.05,  196.19,  195.26,  194.21,  196.31, 194.59, 194.95, 194.26,
      193.06,  194.97,  197.72,  196.07,  195.60,  199.63,  195.30,  191.97,  199.69,  198.34, 196.06, 197.30, 196.43,
      196.07,  195.15,  196.99,  202.76,  203.14,  195.08,  199.69,  202.80,  200.81,  198.33, 203.34, 204.34, 197.69,
      201.61,  197.72,  205.31,  200.97,  200.56,  201.52,  200.08,  207.65,  195.66,  205.47, 202.64, 200.01, 200.90,
      200.07};
   std::vector<Float_t> protonSigma{
      5834.30, 5463.65, 5135.04, 4974.02, 4738.06, 4659.66, 4344.44, 4001.34, 3321.36, 2866.88, 2546.12, 2227.61,
      2018.41, 1900.60, 1822.92, 1706.99, 1562.08, 1441.58, 1409.41, 1369.35, 1304.94, 1192.38, 1082.66, 1039.39,
      1030.27, 1011.92, 967.96,  889.66,  805.44,  779.74,  776.36,  774.98,  767.76,  738.21,  681.22,  624.61,
      601.27,  601.09,  604.94,  605.25,  604.08,  590.53,  559.30,  520.43,  499.05,  495.32,  495.82,  495.89,
      497.52,  494.00,  487.37,  474.17,  451.02,  429.55,  419.49,  415.85,  414.80,  415.23,  412.60,  412.08,
      410.63,  405.80,  396.70,  384.24,  370.44,  360.31,  354.42,  354.02,  353.99,  352.19,  350.52,  348.99,
      346.34,  344.66,  340.97,  335.49,  329.48,  321.13,  313.79,  308.31,  305.99,  305.07,  305.49,  304.03,
      303.40,  301.85,  301.59,  299.92,  298.26,  297.41,  293.96,  291.60,  287.60,  283.08,  278.05,  274.27,
      271.53,  270.70,  268.35,  268.06,  267.85,  267.14,  266.01,  265.72,  264.38,  263.13,  262.19,  260.99,
      259.38,  257.30,  256.34,  253.92,  253.50,  250.85,  251.13,  249.21,  248.20,  245.96,  245.79,  244.22,
      242.20,  240.92,  238.82,  237.26,  236.12,  235.67,  234.11,  232.93,  233.06,  231.74,  231.12,  231.64,
      230.17,  229.69,  229.22,  227.97,  227.48,  227.11,  226.48,  226.15,  225.59,  224.21,  224.02,  222.63,
      222.56,  221.45,  220.12,  220.40,  218.58,  218.86,  217.48,  216.39,  214.97,  215.78,  214.56,  214.06,
      213.02,  213.22,  213.61,  212.78,  212.05,  211.11,  211.98,  209.97,  210.49,  210.29,  210.48,  208.62,
      208.31,  207.24,  207.73,  207.45,  206.01,  205.58,  205.71,  206.14,  205.70,  204.38,  204.30,  204.22,
      204.89,  203.99,  203.72,  203.63,  203.22,  204.82,  202.75,  203.20,  204.08,  203.77,  202.71,  203.79,
      201.84,  201.91,  203.54,  200.67,  202.17,  201.98,  201.37,  202.27,  199.40,  199.81,  200.03,  199.79,
      200.81,  199.39,  200.92,  201.15,  199.94,  199.82,  198.49,  200.69,  198.03,  199.58,  200.80,  197.86,
      198.21,  198.48,  197.28,  198.41,  196.78,  196.24,  198.17,  196.89,  198.92,  195.95,  196.78,  197.02,
      197.05,  195.93,  195.33};

   std::vector<Double_t> retVect; // 0-electron, 1-pion, 2-kaon, 3-proton

   Double_t step      = 0.01;
   Double_t pRange[2] = {0.16, 1.2};

   Int_t    iBin      = 0;
   Double_t dedxMean  = 0.;
   Double_t dedxSigma = 1.;

   // Electron
   if (p >= 0.16 && p <= 1.2) {
      iBin = TMath::FloorNint((p - pRange[0]) / step);

      // std::cout << "e iBin: " << iBin << " size: " << electronMean.size() << std::endl;
      dedxMean  = electronMean.at(iBin);
      dedxSigma = electronSigma.at(iBin);
   } else if (p > 1.2) {
      dedxMean  = 3850.8 + p * 336.062 + p * p * (-115.97);
      dedxSigma = 242.377 + p * 92.0334 + p * p * (-33.1556);
   } else {
      dedxMean  = 3495.38 + p * 2068.7 + p * p * (-2379.39);
      dedxSigma = 233.202 + p * 191.953 + p * p * (-244.544);
   }
   retVect.push_back((dEdx - dedxMean) / dedxSigma);

   pRange[2] = 2.2;

   // Pion
   if (p >= 0.16 && p <= 2.2) {
      iBin = TMath::FloorNint((p - pRange[0]) / step);
      // std::cout << "pi iBin: " << iBin << " size: " << pionMean.size() << std::endl;
      dedxMean  = pionMean.at(iBin);
      dedxSigma = pionSigma.at(iBin);
   } else if (p > 2.2) {
      dedxMean  = 2631.15 + p * 302.052 + p * p * (-36.1923);
      dedxSigma = 174.48 + p * 20.5397 + p * p * (-0.689454);
   } else {
      dedxMean  = 10307.1 + p * (-56934.1) + p * p * 114991;
      dedxSigma = 916.58 + p * (-5205.87) + p * p * 9524.84;
   }
   retVect.push_back((dEdx - dedxMean) / dedxSigma);

   // Kaon
   if (p >= 0.16 && p <= 2.2) {
      iBin = TMath::FloorNint((p - pRange[0]) / step);
      // std::cout << "K iBin: " << iBin << " size: " << kaonMean.size() << std::endl;
      dedxMean  = kaonMean.at(iBin);
      dedxSigma = kaonSigma.at(iBin);
   } else if (p > 2.2) {
      dedxMean  = 3455.91 + p * (-565.856) + p * p * 139.353;
      dedxSigma = 273.621 + p * (-88.8921) + p * p * 25.0721;
   } else {
      dedxMean  = 68436.9 + p * (-382891) + p * p * 627999;
      dedxSigma = 11429.8 + p * (-77924.7) + p * p * 142956;
   }
   retVect.push_back((dEdx - dedxMean) / dedxSigma);

   // Proton
   if (p >= 0.2 && p <= 2.2) {
      iBin = TMath::FloorNint((p - pRange[0]) / step);
      // std::cout << "p iBin: " << iBin << " size: " << protonMean.size() << std::endl;
      dedxMean  = protonMean.at(iBin);
      dedxSigma = protonSigma.at(iBin);
   } else if (p > 2.2) {
      dedxMean  = 4238.02 + p * (-788.916) + p * p * 119.585;
      dedxSigma = 313.054 + p * (-83.9304) + p * p * 14.9712;
   } else {
      dedxMean  = 78147.3 + p * (-258875) + p * p * 253721;
      dedxSigma = 24925.3 + p * (-123173) + p * p * 161350;
   }
   retVect.push_back((dEdx - dedxMean) / dedxSigma);

   return retVect;
}

void MpdMiniDstFillTask::findMatchedMcTracks()
{
   if (!fMCTracks) return;
   fMapRec.resize(fMCTracks->GetEntriesFast());
   std::fill(fMapRec.begin(), fMapRec.end(), false);
   /// mark tracks from E-cal as good
   if (fIsUseECal && fEmcClusters) {
      // Check if MC track that was used in ECal clusters
      for (Int_t iMcTrk = 0; iMcTrk < fMCTracks->GetEntriesFast(); iMcTrk++) {
         for (Int_t iCluster = 0; iCluster < fEmcClusters->GetEntriesFast(); iCluster++) {
            // Retrieve barrel ECal cluster
            MpdEmcClusterKI *cluster = (MpdEmcClusterKI *)fEmcClusters->UncheckedAt(iCluster);

            if (!cluster) continue;

            // Loop over tracks in the cluster
            for (Int_t iTrk = 0; iTrk < cluster->GetNumberOfTracks(); iTrk++) {
               Int_t   id   = -1;
               Float_t eDep = -1.;
               cluster->GetMCTrack(iTrk, id, eDep);
               if (id == iMcTrk) {
                  fMapRec[iMcTrk]            = true;
                  fMcTrk2EcalCluster[iMcTrk] = iCluster;
                  break;
               }
            } // for (Int_t iTrk=0; iTrk < cluster->GetNumberOfTracks(); iTrk++)
         }    // for (Int_t iCluster = 0; iCluster < fEmcClusters->GetEntriesFast(); iCluster++)
      }       // loop over MC tracks
   }

   // store all primary particles
   for (Int_t iMcTrk = 0; iMcTrk < fMCTracks->GetEntriesFast(); iMcTrk++) {
      MpdMCTrack *mcTrack = (MpdMCTrack *)fMCTracks->UncheckedAt(iMcTrk);
      if (mcTrack->GetMotherId() == -1) fMapRec[iMcTrk] = true;
   }

   // store particles matched with kalman, also secondaries
   if (fStoreMatched) {

      for (Int_t iKalmanTrk = 0; iKalmanTrk < fTpcTracks->GetEntriesFast(); iKalmanTrk++) {
         // Retrieve i-th TpcKalmanTrack
         MpdTpcKalmanTrack *kalTrack = (MpdTpcKalmanTrack *)fTpcTracks->UncheckedAt(iKalmanTrk);
         int                id       = kalTrack->GetTrackID();
         fMapRec[id]                 = true;
      }
   }

   // store parents of required MC particles
   if (fStoreParents) {
      for (Int_t iMcTrack = 0; iMcTrack < fMapRec.size(); iMcTrack++) {
         if (fMapRec[iMcTrack]) {
            findMcParent(iMcTrack);
         }
      }
   }
}

void MpdMiniDstFillTask::findMcParent(Int_t index)
{
   MpdMCTrack *track       = (MpdMCTrack *)fMCTracks->UncheckedAt(index);
   Int_t       motherIndex = track->GetMotherId();
   if (motherIndex < 0) return;
   if (motherIndex > fMCTracks->GetEntriesFast()) return;
   if (!fTrivDecay) { // reject "trivial decays:
      MpdMCTrack *mother = (MpdMCTrack *)fMCTracks->UncheckedAt(motherIndex);
      if (mother->GetPdgCode() == 2212) return;                // proton
      if (mother->GetPdgCode() == 2112) return;                // neutron
      if (mother->GetPdgCode() == -2212) return;               // anti-proton
      if (mother->GetPdgCode() == -2112) return;               // anti-neutron
      if (mother->GetPdgCode() == 22) return;                  // photon
      if (mother->GetPdgCode() == track->GetPdgCode()) return; // pdg not changed
   }
   fMapRec[motherIndex] = true;
   findMcParent(motherIndex);
}
