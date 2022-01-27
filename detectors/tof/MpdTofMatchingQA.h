//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_TOF_MATCHING_QA_H
#define __MPD_TOF_MATCHING_QA_H

#include <TList.h>
#include <TString.h>
#include <TVector3.h>
#include <TEfficiency.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#include <iostream>
#include <set>
#include <vector>

#include "MpdTofMatching.h"
#include "MpdTofMatchingData.h"
//------------------------------------------------------------------------------------------------------------------------
class Lsets {
   typedef std::set<size_t> Ts;
   std::vector<Ts>          fData; //!

public:
   Lsets(size_t n) { fData.resize(n); }

   void Insert(size_t index, size_t value);
   void Reset();
   void Print(const char *comment = nullptr, std::ostream &os = std::cout) const;
};
//------------------------------------------------------------------------------------------------------------------------
class TClonesArray;
class MpdTofMatchingQA {
   // QA test histos
   TList fList;

   TH2D *htKFTrack, *htKFTrackPoint, *htTrueMatch, *htMisMatch, *hNMcKfTracks;
   TH2D *hDeltaHitTrue, *hDeltaHitMis, *hDeltaPoint, *hDeltaPoint_Dev, *hDeltaPoint_dR, *hDeltaPoint_dZ,
      *hDeltaPoint_dPhi;
   TH2D *hPvsPt, *hPvsNp, *hPtvsNp, *hPvsP, *hPtvsPt;
   TH1D *hWeightF, *hWeightT, *hNormWeightF, *hNormWeightT;
   TH2D *hPointInsideDetector, *hSmeared_dPhi, *hSmeared_dZ, *hSmeared_dZPhi, *h1, *h2, *h3;

   TEfficiency *pEff1P, *pEff1Eta,
      *pEff1EtaP; // Efficiency = N true matchings / N tpc kf tracks having TOF hit (algorithm efficiency)
   TEfficiency *pCont1P, *pCont1Eta, *pCont1EtaP; // Contamination = N wrong matchings / ( N true matchings + n wrong
                                                  // matchings)(algorithm contamination)
   TEfficiency *pEff2P, *pEff2Eta, *pEff2EtaP;    // Efficiency = N true / N tpc kf tracks (ALICE TDR TOF 2000)
   TEfficiency *pCont2P, *pCont2Eta,
      *pCont2EtaP; // Contamination = N wrong / N single matching to fired strip (ALICE TDR TOF 2000)

   TEfficiency *pEffEtaPIdeal, *pContEtaPIdeal; // ideal matching by MC data
   TEfficiency *pEffKalman, *pEffMatch, *pEffKalmanP, *pEffMatchP, *pEffKalmanPi, *pEffMatchPi, *pEffKalmanK,
      *pEffMatchK;

   TString fFlnm;
   bool    fIsEndcap;

   const char *mangling(const char *name)
   {
      static TString nm;
      nm = fIsEndcap ? "LsETOF_" : "LsTOF_";
      nm += name;
      return nm.Data();
   }
   void Add(TH1 *hist)
   {
      hist->SetDirectory(nullptr);
      fList.Add(hist);
   }
   void Add(TEfficiency *hist)
   {
      hist->SetDirectory(nullptr);
      fList.Add(hist);
   }

public:
   MpdTofMatchingQA(const char *flnm, bool isEndcap);

   void Finish();
   // Fill deviation of estimated point and hit position.
   void FillHitDeviation(const MpdTofMatchingData &);
   // Fill deviation of estimated point and mc point position.
   void FillMcPointDeviation(const TVector3 &mcPosition, const TVector3 &estPointOnPlate, const TVector3 &MCmom,
                             const TVector3 &Momentum, int NofTrHits);
   void FillMatchingEfficiency(const TmPt &tids, const TClonesArray *aTofMatching, const TClonesArray *aTofHits,
                               const TClonesArray *aMCTracks);
   void FillWeightTest(double weight, size_t mult, const TVector3 &, const TVector3 &, const TVector3 &, double angle,
                       const TVector3 &);
   void FillParameter(bool istrue, double weight, double normWeight);
   void FillIdealMatching(const TmPt &mPt, const TmmT2H &mmT2H, const TClonesArray *aMCTracks);
   void FillSmearedPoints(const TVector3 &estPosition, const TVector3 *smearedPosition,
                          const TVector3 *smearedMomentum);
   void FillNtracks(const TClonesArray *aTracks, size_t NkfTracks);

   void FillTrackEfficiency(bool kalman, bool matched, int pdgcode, const TVector3 &momentum);

   inline void FillKfTracks(bool hadPoint, Double_t eta, Double_t Momentum)
   {
      htKFTrack->Fill(eta, Momentum);
      if (hadPoint) htKFTrackPoint->Fill(eta, Momentum);
   }
   inline void FillPointInsideDetector(const TVector3 &position)
   {
      hPointInsideDetector->Fill(position.Z(), position.Phi() * TMath::RadToDeg());
   }

   // Lsets tidsTofTouch = 4;
};
//------------------------------------------------------------------------------------------------------------------------
#endif // # ifndef __MPD_TOF_MATCHING_QA_H
