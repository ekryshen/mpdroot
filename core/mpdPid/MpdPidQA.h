//------------------------------------------------------------------------------------------------------------------------
#ifndef MPD_PID_QA_H
#define MPD_PID_QA_H

//------------------------------------------------------------------------------------------------------------------------
/// \class MpdPidQA
///
/// \brief
/// \author Alexander Mudrokh (LHEP, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------
#include "MpdPid.h"
#include "Rtypes.h"
#include "TGraphAsymmErrors.h"
//------------------------------------------------------------------------------------------------------------------------
using namespace std;
typedef vector<TH1D *> vecTH1Dptrs;

class MpdPidQA : public MpdPid {
public:
   MpdPidQA(); /// default ctor

   MpdPidQA(Double_t sigmaTof, Double_t sigmaEloss,
            Double_t sqrts, /// generators: "URQMD", "LAQGSM" ("QGSM"), "DEFAULT", "NSIG" (for n-sigma method)
            Double_t EnLossCoef = 1., TString Generator = "DEFAULT",
            TString Tracking = "CFHM",    /// tracking: "HP" (Hit Producer), "CF" (Cluster Finder)
            TString NSigPart = "pikapr"); /// possible expressions: el, mu, pi, ka, pr, de, tr, he3, he4

   virtual ~MpdPidQA(); /// destructor

   void                   FillDedxHists(Double_t, Double_t, Int_t);
   void                   Fillm2Hists(Double_t, Double_t, Int_t);
   void                   FillAmplHists(Double_t, Int_t);
   Bool_t                 FillEffContHists(Double_t, Double_t, Int_t, Int_t, Double_t fProbCut = 0.);
   Bool_t                 FillEffContHists(Double_t, Double_t, Double_t, Int_t, Int_t, Double_t fProbCut = 0.);
   void                   GetDedxQA(TString);
   void                   Getm2QA(TString);
   void                   GetAmplQA(TString);
   void                   GetEffContQA(TString);
   MpdPidUtils::ePartType GetPartType(Int_t);

private:
   /// private functions
   void     Init(TString);
   void     FillEffDenominator(Double_t, Int_t);
   Bool_t   FillEffContHists(Double_t, Int_t, Int_t, MpdPidUtils::ePartCharge);
   Double_t Novosibirsk(Double_t *, Double_t *);

   /// variables
   Double_t Xlow[MpdPidUtils::nQAHists];
   Double_t Xhigh[MpdPidUtils::nQAHists];
   Double_t X[MpdPidUtils::nQAHists];
   TString  nSigPart;

   map<Int_t, MpdPidUtils::ePartType> fPartTypeMap; ///< map of correspondence of PDG codes and ePartTypes
   map<MpdPidUtils::ePartType, MpdPidUtils::dEdXStruct_t> fEnLossMap;          ///< dE/dx QA map
   map<MpdPidUtils::ePartType, TH2D *>                    fMSquaredMap;        ///< m^2 QA map
   map<MpdPidUtils::ePartType, vecTH1Dptrs>               fAbundanceMap;       ///< Abundance QA map
   map<MpdPidUtils::ePartType, std::vector<vecTH1Dptrs>>  fEffContMap;         ///< PID efficiency and contamination map
   TH2D                                                  *fSumBetheBlochHists; ///< dE/dx VS. p histogram, all species
   TH2D *fChBetheBlochHists; ///< dE/dx VS. p histogram, neg char species plotted with p -> (-p)
   TH2D *fm2LightHist;       ///< m^2 VS. p histogram for pi, K and (anti-)p
   TH2D *fm2HeavyHist;       ///< m^2 VS. p histogram for d, t and he-3

   ClassDef(MpdPidQA, 3);
};

#endif
