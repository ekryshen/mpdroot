/**
 * \class MpdMiniDstFillTrask
 * \brief Class that reads event, track, and detector parameters and fills MpdMiniDst
 *
 * The MpdMiniDstFillTask reads information from the reconstructed (and/or simulated) event, track, vertex
 * and detectors and fill MpdMiniDst structure.
 *
 * \author Pavel Batyuk (JINR), Grigory Nigmatkulov (NRNU MEPhI)
 * \email pavel.batyuk@jinr.ru ; nigmatkulov@gmail.com ; ganigmatkulov@mephi.ru
 * \date Aug 01, 2019
 **/

#ifndef MpdMiniDstFillTask_h
#define MpdMiniDstFillTask_h

// C++ headers
#include <vector>
#include <utility>
#include <map>

// ROOT and FAIRROOT headers
#include "TNamed.h"
#include "TClass.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TClonesArray.h"
#include "TFile.h"

// FairRoot headers
#include "FairTask.h"
#include "FairRunAna.h"
#include "FairField.h"
#include "FairRootManager.h"
#include "FairMCEventHeader.h"
#include "FairEventHeader.h"

// MpdRoot headers
class MpdEvent;
class MpdKalmanTrack;
class MpdVertex;
class MpdTpcKalmanTrack;

// MpdMiniDst headers
class MpdMiniDst;
class MpdMiniTrack;
class MpdMiniTrackCovMatrix;

//_________________
class MpdMiniDstFillTask : public FairTask {

public:
   /// Default constructor
   MpdMiniDstFillTask();

   /// Constructor with output file name for miniDST and IO-modes (read or write)
   MpdMiniDstFillTask(TString);
   /// Destructor
   virtual ~MpdMiniDstFillTask();

   /// Init
   virtual InitStatus Init();
   /// Execute
   virtual void Exec(Option_t *option);
   /// Finish
   virtual void Finish();

   /// Fill track covariant matrix
   /// \param true store covariant matrix (default)
   /// \param false do not store covariant matrix information
   void isUseCovMatrix(Bool_t flag) { fIsUseCovMatrix = flag; }

   /// Fill ECal-related information
   void isUseECal(Bool_t flag) { fIsUseECal = flag; }
   /// Fill all MC particles that are reconstructed
   void isUseMatched(Bool_t flag) { fStoreMatched = flag; }
   /// Fill all MC particles that are parent of other needed MC particles
   void isUseParent(Bool_t flag) { fStoreParents = flag; }
   /// Fill/not fill parents from trivial processes e.g. p->n n->p, x->x gamma->x work only with useParent
   void isStoreTrivialDecays(Bool_t flag) { fTrivDecay = flag; };

   ///  Set how to fill nSigma(e,pi,K,p)
   /// \par 0 use values from MpdDst (default)
   /// \par 1 use values from special tables (cheat way)
   void setNSigmaDedxEstimator(Int_t estimator)
   {
      fNSigmaDedxEstimator = (estimator == 0 || estimator == 1) ? estimator : 0;
   }

private:
   /// Store/Not store track covariant matrix information
   Bool_t fIsUseCovMatrix;

   /// Store/Not store ECal information
   Bool_t fIsUseECal;
   /// Store/Not secondary mc tracks that were reconstructed
   Bool_t fStoreMatched;
   /// store parents of written particles
   Bool_t fStoreParents;
   // store parents that are stable particles
   Bool_t fTrivDecay;
   /// Pointer to event header
   FairMCEventHeader *fEventHeaders;
   /// Pointer to MpdEvent
   MpdEvent *fEvents;
   /// Pointer to primary vertieces
   TClonesArray *fVertices;
   /// Pointer to TPC tracks
   TClonesArray *fTpcTracks;
   /// Pointer to TOF hits
   TClonesArray *fTofHits;
   /// Pointer to TOF-matching information
   TClonesArray *fTofMatching;
   /// Pointer to MC tracks
   TClonesArray *fMCTracks;
   // Pointer to GenTracks (Primaries from MC-generator)
   TClonesArray *fGenTracks;
   // Pointer to EMC clusters
   TClonesArray *fEmcClusters;
   // Pointer to ZDC digits
   TClonesArray *fZdcDigits;

   /// A pointer to the main input/outpur miniDst structure containing all `TObjArray`s
   MpdMiniDst *fMiniDst;

   /// Magnetic field of the current event
   Float_t fBField;

   /// Output file name
   TString fOutputFileName;
   /// Pointer to the output file
   TFile *fOutputFile;

   /// Pointer to the TTree with miniDst
   TTree *fTTree;

   /// Splitting level of ROOT file
   Int_t fSplit;
   /// Compression level
   Int_t fCompression;
   /// Size of the buffer
   Int_t fBufferSize;

   /// Pointer to the mini arrays
   TClonesArray **fMiniArrays;

   /// Vector that keeps McTrack to miniMcTrack correspondence
   std::vector<std::pair<Int_t, UShort_t>> fMcTrk2MiniMcTrk;
   /// Mat for keeping MC track to barrel ECal cluster correpsondence
   std::map<Int_t, Int_t> fMcTrk2EcalCluster;
   // map of reconstructed mc tracks
   std::vector<bool> fMapRec;

   /// How to fill nSigma(e,pi,K,p)
   /// \par 0 use values from MpdDst (default)
   /// \par 1 use values from special tables (cheat way)
   Int_t fNSigmaDedxEstimator;

   /// Turn-off ROOT streamers
   void streamerOff();

   /// Create arrays
   void createArrays();

   /// Fill event information
   void fillEvent();

   /// Fill MC tracks
   void fillMcTracks();

   /// Fill barrel ECal clusters
   void fillECalClusters();

   /// Fill track information
   void fillTracks();

   /// Fill BTOF information
   void fillBTofHits();

   /// Fill FHCal information
   void fillFHCalHits();

   // marks matched mc tracks
   void findMatchedMcTracks();

   /// Return index of miniMcTrack that corresponds to McTrack (-1 not found)
   Int_t miniMcIdxFromMcIdx(Int_t mcIdx);

   /* Below is given a set of auxiliary functions */
   /// Compute matrices of derivatives
   void computeAandB(TMatrixD &, const MpdKalmanTrack *, const MpdKalmanTrack &, TMatrixD &, TMatrixD &, TMatrixD &);

   /// Adjust track parameters
   void Proxim(const MpdKalmanTrack &, MpdKalmanTrack &);

   /// A function to be used when doing refit of tracks considered
   /// as primaries to the primary vertex aimed at precising of their momenta
   void refit2Vp(MpdMiniTrack *, Int_t, MpdVertex *);

   /// Fill covariance matrix if needed
   void fillCovMatrix(MpdTpcKalmanTrack *, MpdMiniTrackCovMatrix *);

   /// find parents (and grandparents of mc particle with given index and mark them in fMapRec
   void findMcParent(Int_t index);

   /// Check if the event event is okay (reasonabler reconstruction)
   Bool_t isGoodEvent();

   /// Estimate nSigma value (cheat way). Returns vector with
   /// values for electron, pion, kaon and proton
   std::vector<Double_t> nSigmaDedx(Double_t p, Double_t dEdx);

   ClassDef(MpdMiniDstFillTask, 0)
};

#endif // MpdMiniDstFillTask_h
