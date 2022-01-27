//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_TOF_T0_H
#define __MPD_TOF_T0_H

//------------------------------------------------------------------------------------------------------------------------
/// \class MpdTofT0
///
/// \brief
/// \author Viktor Baryshnikov (Faculty of Physics, MSU, Moscow)
/// \author Sergei Lobastov (LHE, JINR, Dubna)
// usage:
// insert MpdTofT0 task into FairRunAna task chain AFTER MpdTofMatching task.
//	auto tofT0 = new MpdTofT0("task name", Bool_t useMCdata = true, Int_t verbose = 1, "QAflnm.root")
// or 	auto tofT0 = new MpdTofT0; // default ctr for MC data run and without QA tests
//	fRun->AddTask(tofT0) -- add task to chain
// result:
// Tevent - T event by TOF [ns],  Chi2 - weight, nMatchings - number of samples used
// public members of MpdTofT0Data class.
//  FairRootManager::Instance()->Register("TOFTevent", "Tof", MpdTofT0Data pointer, kFALSE);
//------------------------------------------------------------------------------------------------------------------------
#include <vector>
#include <map>

#include <TList.h>
#include <TString.h>
#include <TNamed.h>

#include "FairTask.h"
//------------------------------------------------------------------------------------------------------------------------
class TRandom2;
class TClonesArray;
class TStopwatch;
class TH1D;
class TH2D;
class TEfficiency;

class FairMCEventHeader;
//------------------------------------------------------------------------------------------------------------------------
class LState {
   const static size_t fBase = 3; // pion, proton, kaon
   size_t              fSize;
   const double        fMass[fBase]   = {0.13957, 0.938272, 0.493677}; // [GeV] species mass
   const double        fTsigma[fBase] = {
      10., 10., 10.}; // [ns]  uncertainty (sigma_Texp,i) due to the tracking and reconstruction MUST BE CORRECTED!!!
   std::vector<char> fState; // char = [0,1,2] == [pion, proton, kaon]

   void _checkAndShift(size_t index);

public:
   LState(size_t &size); // max value = 19 matching, pow(3,19) = 1,162,261,500 combination of hypothesis(~10Gb)
   ~LState();

   char  &operator[](size_t index) { return fState[index]; };
   double Mass(size_t index) { return fMass[fState[index]]; };
   double TSigma(size_t index) { return fTsigma[fState[index]]; };

   void   Next();       // switch state to next hypothesis
   size_t Size() const; // number of combination(pow(fBase,size))
   void   Print(const char *comment = nullptr, std::ostream &os = std::cout) const; // print current hypothesis
};
//------------------------------------------------------------------------------------------------------------------------
class MpdTofT0Data : public TNamed {
public:
   double Chi2       = {}; // chi2 or weight  = 1/(Chi2/NDF)
   double Tevent     = {}; // Tevent Tof estimation
   size_t nMatchings = {};

   MpdTofT0Data() : TNamed("", ""){};
   MpdTofT0Data(const char *name) : TNamed(name, ""){};

   void reset()
   {
      Chi2       = 0.;
      Tevent     = 0.;
      nMatchings = 0;
   }

   ClassDef(MpdTofT0Data, 1);
};
//------------------------------------------------------------------------------------------------------------------------
class MpdTofT0 : public FairTask {
// number of matching inside one Pt division.
#define fRangeSize 15

   typedef std::multimap<double, std::pair<MpdTofMatchingData *, MpdTpcKalmanTrack *>> T_matchData; // key = Pt

   typedef struct TmatchMini {
      double                P;    // Momentum [GeV/c]
      double                L;    // Track length [cm]
      double                Ttof; // Tof time [ns]
      double                Texp;
      double                Weight;
      T_matchData::iterator it; // link to data

      TmatchMini() : P(0.), L(0.), Ttof(0.){};
      TmatchMini(double p, double l, double t, T_matchData::iterator i) : P(p), L(l), Ttof(t), it(i){};
   } T_matchMini;

   MpdTofT0Data      *fResult = nullptr;
   T_matchData        mmMatchData;
   double             fTofSigma     = 60.;     //[ns]
   TClonesArray      *aTofMatchings = nullptr; //! <---  input
   TClonesArray      *aTPCkfTracks  = nullptr; //! <---  input
   TClonesArray      *aMcTracks     = nullptr; //! <--- MC input
   FairMCEventHeader *fEventHeader  = nullptr; //! <--- MC input

   TRandom2   *pRandom = nullptr;
   TStopwatch *pTimer  = nullptr;
   TList       fList;
   TString     fFlnm;
   Bool_t      fUseMCData;
   Bool_t      fDoTest = false;

   // QA histos
   TH1D *hNmatch = nullptr;
   TH2D *hTeveN = nullptr, *hTcalcN = nullptr, *hTofmultB = nullptr;
   // efficiency of T event better 10, 15, 20, 50, 100 ps
   TEfficiency *effT0_10 = nullptr, *effT0_15 = nullptr, *effT0_20 = nullptr, *effT0_50 = nullptr, *effT0_100 = nullptr;

   MpdTofT0Data Estimate(T_matchMini *data, size_t size);

public:
   MpdTofT0(const char *name = "TOF_Tevent", Bool_t useMCdata = true, Int_t verbose = 1, const char *flnm = nullptr);
   virtual ~MpdTofT0();

   virtual InitStatus Init();
   virtual void       Exec(Option_t *option);
   virtual void       Finish();

   void SetTofSigma(Double_t sigma) { fTofSigma = sigma; };
   void SetSeed(UInt_t seed = 0);

   template <typename T>
   void Add(T *hist);

   ClassDef(MpdTofT0, 2);
};
//------------------------------------------------------------------------------------------------------------------------
#endif // #ifndef __MPD_TOF_T0_H
