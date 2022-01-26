//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_ETOF_MATCHING_H
#define __MPD_ETOF_MATCHING_H 1

//------------------------------------------------------------------------------------------------------------------------
/// \class MpdEtofMatching
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------

#include <TObject.h>
#include <TVector3.h>

#include "MpdTofUtils.h"
#include "MpdTofHit.h"
#include "MpdTofMatchingData.h"

#include "FairTask.h"

//------------------------------------------------------------------------------------------------------------------------
class MpdEctKalmanTrack;
class TClonesArray;
// class LMatchingFilter;
class MpdTofMatchingQA;
class MpdKalmanFilter;
//------------------------------------------------------------------------------------------------------------------------
class MpdEtofMatching : public FairTask {
public:
   enum kMode { kBruteForce, kIntervalTree, kSmearingDummyTracks };

private:
   TClonesArray *aMcPoints;     //! <--- MC input
   TClonesArray *aMcTracks;     //! <--- MC input
   TClonesArray *aTofHits;      //! <--- input TOF hits
   TClonesArray *aKFectTracks;  //! <--- input KF Ect Tracks
   TClonesArray *aTofMatchings; //! ---> output

   kMode     fMode;
   Bool_t    fDoTest, fUseMCData, fDoMCTest;
   TRandom2 *pRandom;
   Int_t     fNSmeared;              // default = 20, N samples for kSmearingDummyTracks method
   Double_t  fTofZpos, fTofRmax;     // [cm]
   Double_t  fThreshR, fThreshTheta; // [cm], acceptance window

   MpdKalmanFilter * pKF;         //!
   MpdTofMatchingQA *pMatchingQA; //!
                                  //	LMatchingFilter		*pMF; 		//!

   bool EstTrackOnPlane(const MpdEctKalmanTrack *tr, Double_t Z, TVector3 &pos, Double_t &length, TVector3 &Mom,
                        Int_t &charge) const;

public:
   MpdEtofMatching(const char *name = "ETOF Matching", Int_t verbose = 1, Bool_t DoTest = false,
                   const char *flnm = "QA.MpdEtofMatching.root");
   ~MpdEtofMatching();

   InitStatus Init();
   void       Exec(Option_t *option);
   void       Finish();

   void SetETofZposition(Double_t z) { fTofZpos = z; };

   ClassDef(MpdEtofMatching, 1) // MpdEtofMatching
};
//------------------------------------------------------------------------------------------------------------------------
#endif
