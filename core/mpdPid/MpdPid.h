//------------------------------------------------------------------------------------------------------------------------
#ifndef MPD_PID_H
#define MPD_PID_H

//------------------------------------------------------------------------------------------------------------------------
/// \class MpdPid
///
/// \brief
/// \author Alexander Mudrokh (LHEP, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------
#include "MpdPidUtils.h"
#include "MpdTrack.h"
//#include "MpdTpcKalmanTrack.h"  /// MiniDST
#include "FairRunAna.h" /// MiniDST
#include <TFile.h>
#include <TF2.h>
#include <TMatrixD.h> /// MiniDST
//------------------------------------------------------------------------------------------------------------------------
using namespace std;
typedef vector<TF1 *> vecTF1ptrs;

class MpdPid : public TObject {
public:
   //--------------------------------------------------------------------//
   //--------------------- Constructor / Destructor ---------------------//
   //--------------------------------------------------------------------//

   MpdPid(); /// default ctor

   MpdPid(Double_t sigmaTof, Double_t sigmaEloss,
          Double_t sqrts, /// generators: "PHSD", "URQMD", "LAQGSM" ("QGSM"), "DEFAULT", "EPOS" (for pp collisions),
                          /// "NSIG" (for n-sigma method)
          Double_t EnLossCoef = 1., TString Generator = "DEFAULT",
          TString Tracking = "CFHM",    /// tracking: "HP" (Hit Producer), "CF" (Cluster Finder)
          TString NSigPart = "pikapr"); /// possible expressions: pi, ka, pr, el, mu, de, tr, he3, he4

   virtual ~MpdPid(); /// destructor

   //--------------------------------------------------------------------//
   //---------------------- Array of probabilities ----------------------//
   //--------------------------------------------------------------------//

   /// Fill array of probabilities, otherwise return kFALSE
   Bool_t FillProbs(MpdTrack *);
   Bool_t FillProbs(Double_t, Double_t, Int_t);           /// variables: full momentum, dE/dx, charge
   Bool_t FillProbs(Double_t, Double_t, Double_t, Int_t); /// variables: full momentum, dE/dx, mass squared, charge

   /// Return probabilities
   Double_t GetProbEl(void) { return fProb[MpdPidUtils::kElectron]; }
   Double_t GetProbMu(void) { return fProb[MpdPidUtils::kMuon]; }
   Double_t GetProbPi(void) { return fProb[MpdPidUtils::kPion]; }
   Double_t GetProbKa(void) { return fProb[MpdPidUtils::kKaon]; }
   Double_t GetProbPr(void) { return fProb[MpdPidUtils::kProton]; }
   Double_t GetProbDe(void) { return fProb[MpdPidUtils::kDeuteron]; }
   Double_t GetProbTr(void) { return fProb[MpdPidUtils::kTriton]; }
   Double_t GetProbHe3(void) { return fProb[MpdPidUtils::kHe3]; }
   Double_t GetProbHe4(void) { return fProb[MpdPidUtils::kHe4]; }
   Double_t GetProb(MpdPidUtils::ePartType iType) { return fProb[iType]; }

   Long_t GetMaxProb(); /// Returns the most probable pdg code

   //--------------------------------------------------------------------//
   //------------------------------ dE/dx -------------------------------//
   //--------------------------------------------------------------------//

   /// Returns expected <dE/dx> value (variable: full momentum)
   Double_t GetDedxElParam(Double_t p) { return GetDedxParam(p, MpdPidUtils::kElectron); }
   Double_t GetDedxMuParam(Double_t p) { return GetDedxParam(p, MpdPidUtils::kMuon); }
   Double_t GetDedxPiParam(Double_t p) { return GetDedxParam(p, MpdPidUtils::kPion); }
   Double_t GetDedxKaParam(Double_t p) { return GetDedxParam(p, MpdPidUtils::kKaon); }
   Double_t GetDedxPrParam(Double_t p) { return GetDedxParam(p, MpdPidUtils::kProton); }
   Double_t GetDedxDeParam(Double_t p) { return GetDedxParam(p, MpdPidUtils::kDeuteron); }
   Double_t GetDedxTrParam(Double_t p) { return GetDedxParam(p, MpdPidUtils::kTriton); }
   Double_t GetDedxHe3Param(Double_t p) { return GetDedxParam(p, MpdPidUtils::kHe3); }
   Double_t GetDedxHe4Param(Double_t p) { return GetDedxParam(p, MpdPidUtils::kHe4); }
   Double_t GetDedxParam(Double_t, MpdPidUtils::ePartType);

   /// Returns expected dE/dx width (variables: full momentum, species)
   Double_t GetDedxWidthValue(Double_t, MpdPidUtils::ePartType);

   /// Returns expected asymmetry value delta (variables: full momentum, species)
   Double_t GetTailValue(Double_t, MpdPidUtils::ePartType);

   Double_t GetdEdxBottomBound(MpdPidUtils::ePartType iType) { return fdEdxBottomBound[iType]; }
   Double_t GetdEdxTopBound(MpdPidUtils::ePartType iType) { return fdEdxTopBound[iType]; }
   void     SetdEdxBottomBound(Double_t nsig, MpdPidUtils::ePartType iType) { fdEdxBottomBound[iType] = nsig; }
   void     SetdEdxTopBound(Double_t nsig, MpdPidUtils::ePartType iType) { fdEdxTopBound[iType] = nsig; }

   /// Returns expexted <dE/dx>, dE/dx width and asymmetry parameterizations as vectors of TF1*
   vecTF1ptrs GetVecdEdxMean(MpdPidUtils::ePartType);
   vecTF1ptrs GetVecdEdxWidth(MpdPidUtils::ePartType);
   vecTF1ptrs GetVecdEdxAsym(MpdPidUtils::ePartType);

   //--------------------------------------------------------------------//
   //------------------------------- m^2 --------------------------------//
   //--------------------------------------------------------------------//

   /// Returns expected m^2 width (variables: full momentum, species)
   Double_t Getm2WidthParam(Double_t, MpdPidUtils::ePartType);

   /// Returns expected m^2 width parameterization as vector of TF1*
   vecTF1ptrs GetVecm2Width(MpdPidUtils::ePartType);

   Double_t Getm2BottomBound(MpdPidUtils::ePartType iType) { return fm2BottomBound[iType]; }
   Double_t Getm2TopBound(MpdPidUtils::ePartType iType) { return fm2TopBound[iType]; }
   void     Setm2BottomBound(Double_t nsig, MpdPidUtils::ePartType iType) { fm2BottomBound[iType] = nsig; }
   void     Setm2TopBound(Double_t nsig, MpdPidUtils::ePartType iType) { fm2TopBound[iType] = nsig; }

   //--------------------------------------------------------------------//
   //----------------------------- Yileds -------------------------------//
   //--------------------------------------------------------------------//

   Double_t GetBayesCoefficient(Double_t, MpdPidUtils::ePartType, MpdPidUtils::ePartCharge);

   Double_t GetPrRat() { return fPrRatio; } /// Returns pos./neg. ratio for protons
   void     SetPrRat(Double_t);             /// Set pos./neg. ratio for protons

   /// Returns expected yield parameterization as vector of TF1*
   /// use via vec[MpdPidUtils::kPos/kNeg]
   vecTF1ptrs GetVecYield(MpdPidUtils::ePartType);

   //--------------------------------------------------------------------//
   //------------------------------ Other -------------------------------//
   //--------------------------------------------------------------------//

   void Print(const char *comment, ostream &os);

   /// Returns distance to the most probable dE/dx and m2 value in terms of sigmas
   /// possible expressions: pi, ka, pr, el, mu, de, tr, he3, he4
   /// use it after MpdPid::FillProbs usage only!!!
   Double_t GetNsigmaToBetheBloch(TString);
   Double_t GetNsigmaToBetheBloch(MpdPidUtils::ePartType);
   Double_t GetNsigmaToAverageMass2(TString);
   Double_t GetNsigmaToAverageMass2(MpdPidUtils::ePartType);

   /// Includes particle species iType to n-sigma method
   void SetParticleEnabledNsig(MpdPidUtils::ePartType iType) { fNSigSpecies[iType] = kTRUE; }

   /// Additional functions
   Double_t parElBB(Double_t *, Double_t *);
   Double_t parMuBB(Double_t *, Double_t *);
   Double_t parPiBB(Double_t *, Double_t *);
   Double_t parKaBB(Double_t *, Double_t *);
   Double_t parPrBB(Double_t *, Double_t *);
   Double_t parDeBB(Double_t *, Double_t *);
   Double_t parTrBB(Double_t *, Double_t *);
   Double_t parHe3BB(Double_t *, Double_t *);
   Double_t parHe4BB(Double_t *, Double_t *);

   Double_t AsymGaus(Double_t *, Double_t *);
   Double_t AsymGaus2(Double_t *, Double_t *);
   Double_t MomPi(Double_t *, Double_t *);
   Double_t MomPr(Double_t *, Double_t *);

protected:
   void Init(TString, TString, TString, Double_t);

   /// GetDedxProb_asym variables: cut_dedx, p, dedx, n, emean, sige, species
   Double_t ComputeDedxProb_asym(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, MpdPidUtils::ePartType);

   /// GetCombProb_asym variables: cut_dedx, cut_m2, p, dedx, m2, n, emean, mmean, sige, sigm, species
   Double_t ComputeCombProb_asym(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t,
                                 Double_t, Double_t, MpdPidUtils::ePartType);

   void ComputeBayesCoefficients(Double_t);

   Double_t ComputeEnLossSigma(Double_t, Double_t, Double_t, Double_t, MpdPidUtils::ePartType);
   Double_t ComputeMSquaredSigma(Double_t, Double_t, Double_t, MpdPidUtils::ePartType);

   /// Variables

   TF1 *fGaus;
   TF1 *fAsymGaus;
   TF2 *fGaus2;
   TF2 *fAsymGaus2;

   Double_t fProb[MpdPidUtils::kNSpecies + 1]; ///< the probability to identify track as a species <i>
   Double_t
      fEnLossSigmasArray[MpdPidUtils::kNSpecies]; ///< the deviation of the measured energy loss from that expected for
                                                  ///< the species <i>, in terms of the detector resolution
   Double_t
      fMSquaredSigmasArray[MpdPidUtils::kNSpecies];   ///< the deviation of the measured mass squared from that expected
                                                      ///< for the species <i>, in terms of the detector resolution
   Double_t fdEdxBottomBound[MpdPidUtils::kNSpecies]; ///< dE/dx bottom bound ( dE/dx < <dE/dx> ) in terms of sigmas
   Double_t fdEdxTopBound[MpdPidUtils::kNSpecies];    ///< dE/dx top bound ( dE/dx > <dE/dx> ) in terms of sigmas
   Double_t fm2BottomBound[MpdPidUtils::kNSpecies];   ///< m^2 bottom bound ( dE/dx < <dE/dx> ) in terms of sigmas
   Double_t fm2TopBound[MpdPidUtils::kNSpecies];      ///< m^2 top bound ( dE/dx > <dE/dx> ) in terms of sigmas
   Double_t fPrRatio;                                 ///< proton ratio is pos./neg.
   Double_t fEnergy;                                  ///< collision energy
   MpdPidUtils::eTrackingState fTrackingState;
   Bool_t                      fMethod; ///< PID method, kTRUE - Bayesian approach, kFALSE - n-sigma method
   MpdPidUtils::ePartCharge    fCharge; ///< track charge
   Double_t                    fBayesCoefficients[MpdPidUtils::kNSpecies]; ///< Bayes Coefficients assigned to the track
   Bool_t fNSigSpecies[MpdPidUtils::kNSpecies]; ///< array of flags to the particle species involved in n-sigma method

   map<MpdPidUtils::ePartType, vecTF1ptrs> fdEdxBBMap;    /// Bethe-Bloch functions for mean energy deposit description
   map<MpdPidUtils::ePartType, vecTF1ptrs> fdEdxSigmaMap; /// Polynomials for dE/dx width description
   map<MpdPidUtils::ePartType, vecTF1ptrs> fdEdxDeltaMap; /// Polynomials for dE/dx asymmetry description
   map<MpdPidUtils::ePartType, vecTF1ptrs> fParM2Map;     /// Functions for m2 width description
   map<MpdPidUtils::ePartType, vecTF1ptrs> fPartYieldMap; /// Functions for multiplicity description

   ClassDef(MpdPid, 4);
};

#endif
