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

   /*!
     \param sigmaTof
     \param sigmaEloss
     \param sqrts
     \param EnLossCoef
     \param Generator <table>
   <tr><td>"URQMD"</td></tr>
   <tr><td>"LAQGSM" ("QGSM")</td></tr>
   <tr><td>"DEFAULT"</td></tr>
   <tr><td>"NSIG" (for n-sigma method)</td></tr>
   </table>
     \param Tracking <table>
   <tr><td>"HP" (Hit Producer)</td></tr>
   <tr><td>"CF" (Cluster Finder)</td></tr>
   </table>
     \param NSigPart <table>
   <tr><td>pi</td></tr>
   <tr><td>ka</td></tr>
   <tr><td>pr</td></tr>
   <tr><td>el</td></tr>
   <tr><td>mu</td></tr>
   <tr><td>de</td></tr>
   <tr><td>tr</td></tr>
   <tr><td>he3</td></tr>
   <tr><td>he4</td></tr>
   </table>
   */
   MpdPidQA(Double_t sigmaTof, Double_t sigmaEloss, Double_t sqrts, Double_t EnLossCoef = 1.,
            TString Generator = "DEFAULT", TString Tracking = "CFHM", TString NSigPart = "pikapr");

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
