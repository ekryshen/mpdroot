#include "MpdPidQA.h"

MpdPidQA::MpdPidQA() : MpdPid()
{
   /// Default constructor
   fSumBetheBlochHists = nullptr;
   fChBetheBlochHists  = nullptr;
   fm2LightHist        = nullptr;
   fm2HeavyHist        = nullptr;
   Init("pikapr");
}

MpdPidQA::MpdPidQA(Double_t sigmaTof, Double_t sigmaEloss, Double_t sqrts, Double_t EnLossCoef, TString Generator,
                   TString Tracking, TString NSigPart)
   : MpdPid(sigmaTof, sigmaEloss, sqrts, EnLossCoef, Generator, Tracking, NSigPart)
{
   Init(NSigPart);
}

MpdPidQA::~MpdPidQA()
{
   for (auto it = fEnLossMap.begin(); it != fEnLossMap.end(); ++it) {
      for (Int_t i = 0; i < MpdPidUtils::nQAHists; i++) delete it->second.dEdXPart[i];
      delete it->second.BetheBlochHist;
   }
   for (auto it = fMSquaredMap.begin(); it != fMSquaredMap.end(); ++it) delete it->second;
   for (auto it = fAbundanceMap.begin(); it != fAbundanceMap.end(); ++it) {
      for (TH1D *hYield : it->second) {
         delete hYield;
      }
   }
   for (auto it = fEffContMap.begin(); it != fEffContMap.end(); ++it) {
      for (vecTH1Dptrs vecEC : it->second) {
         for (TH1D *hEC : vecEC) {
            delete hEC;
         }
      }
   }
   delete fSumBetheBlochHists;
   delete fChBetheBlochHists;
   delete fm2LightHist;
   delete fm2HeavyHist;
}

void MpdPidQA::Init(TString Particles)
{
   cout << "MpdPidQA::Init().." << endl;
   Double_t    PMIN = 0.0, PMAX = 5.0;
   Double_t    XMIN[MpdPidUtils::nQAHists][MpdPidUtils::kNSpecies], XMAX[MpdPidUtils::nQAHists][MpdPidUtils::kNSpecies];
   TString     histName;
   const Int_t nbins = 200, nbins2 = 100;
   for (Int_t i = 0; i < MpdPidUtils::nQAHists; i++) {
      for (Int_t j = 0; j < MpdPidUtils::kNSpecies; j++) {
         XMIN[i][j] = -1.;
         XMAX[i][j] = -1.;
      }
   }
   TString ECHistNames[4] = {"All", "Id", "IdRight", "IdWrong"};
   nSigPart               = Particles;
   Double_t MinEnLoss = 0.0, MaxEnLoss = 0.0, MaxEnLossLoc;

   /// filling out fPartTypeMap
   fPartTypeMap.insert(pair<Int_t, MpdPidUtils::ePartType>(11, MpdPidUtils::kElectron));
   fPartTypeMap.insert(pair<Int_t, MpdPidUtils::ePartType>(-11, MpdPidUtils::kElectron));
   fPartTypeMap.insert(pair<Int_t, MpdPidUtils::ePartType>(13, MpdPidUtils::kMuon));
   fPartTypeMap.insert(pair<Int_t, MpdPidUtils::ePartType>(-13, MpdPidUtils::kMuon));
   fPartTypeMap.insert(pair<Int_t, MpdPidUtils::ePartType>(211, MpdPidUtils::kPion));
   fPartTypeMap.insert(pair<Int_t, MpdPidUtils::ePartType>(-211, MpdPidUtils::kPion));
   fPartTypeMap.insert(pair<Int_t, MpdPidUtils::ePartType>(321, MpdPidUtils::kKaon));
   fPartTypeMap.insert(pair<Int_t, MpdPidUtils::ePartType>(-321, MpdPidUtils::kKaon));
   fPartTypeMap.insert(pair<Int_t, MpdPidUtils::ePartType>(2212, MpdPidUtils::kProton));
   fPartTypeMap.insert(pair<Int_t, MpdPidUtils::ePartType>(-2212, MpdPidUtils::kProton));
   fPartTypeMap.insert(pair<Int_t, MpdPidUtils::ePartType>(PDG_DEUTERON, MpdPidUtils::kDeuteron));
   fPartTypeMap.insert(pair<Int_t, MpdPidUtils::ePartType>(PDG_TRITON, MpdPidUtils::kTriton));
   fPartTypeMap.insert(pair<Int_t, MpdPidUtils::ePartType>(PDG_HE3, MpdPidUtils::kHe3));
   fPartTypeMap.insert(pair<Int_t, MpdPidUtils::ePartType>(PDG_HE4, MpdPidUtils::kHe4));

   if (Particles.Contains("el")) {
      /// dE/dx
      if ((fTrackingState == MpdPidUtils::kCF) || (fTrackingState == MpdPidUtils::kCFHM)) {
         MaxEnLossLoc = 5500.;
         for (Int_t i = 0; i < 5; i++) {
            XMIN[i][MpdPidUtils::kElectron] = 2000.;
            XMAX[i][MpdPidUtils::kElectron] = 5000.;
         }
         for (Int_t i = 5; i < MpdPidUtils::nQAHists; i++) {
            XMIN[i][MpdPidUtils::kElectron] = 2500.;
            XMAX[i][MpdPidUtils::kElectron] = 5500.;
         }
      } else {
         MaxEnLossLoc = 3.e-06;
         for (Int_t i = 0; i < MpdPidUtils::nQAHists; i++) {
            XMIN[i][MpdPidUtils::kElectron] = 1.e-06;
            XMAX[i][MpdPidUtils::kElectron] = 3.e-06;
         }
      }
      MpdPidUtils::dEdXStruct_t SEl;
      for (Int_t iHist = 0; iHist < MpdPidUtils::nQAHists; iHist++) {
         histName = "El_" + TString::Itoa(iHist, 10);
         SEl.dEdXPart[iHist] =
            new TH1D(histName, "", nbins, XMIN[iHist][MpdPidUtils::kElectron], XMAX[iHist][MpdPidUtils::kElectron]);
      }
      SEl.BetheBlochHist = new TH2D("hElBB_QA", "", 300, 0.05, 3.0, 300, MinEnLoss, MaxEnLossLoc);
      SEl.ibeg           = 0;
      SEl.iend           = 39;
      fEnLossMap.insert(pair<MpdPidUtils::ePartType, MpdPidUtils::dEdXStruct_t>(MpdPidUtils::kElectron, SEl));

      /// m-squared
      TH2D *mass2El = new TH2D("mass2El", "", nbins, 0., 4., nbins, -0.7, 0.7);
      fMSquaredMap.insert(pair<MpdPidUtils::ePartType, TH2D *>(MpdPidUtils::kElectron, mass2El));

      /// Abundance
      vecTH1Dptrs vecYield;
      TH1D       *hYield_elpos = new TH1D("hYield_elpos", "", 600, 0., 3.);
      vecYield.push_back(hYield_elpos);
      TH1D *hYield_elneg = new TH1D("hYield_elneg", "", 600, 0., 3.);
      vecYield.push_back(hYield_elneg);
      fAbundanceMap.insert(pair<MpdPidUtils::ePartType, vecTH1Dptrs>(MpdPidUtils::kElectron, vecYield));

      /// Efficiency & Contamination
      TH1D                    *hEC_pos[4], *hEC_neg[4];
      vecTH1Dptrs              vecEC_pos, vecEC_neg;
      std::vector<vecTH1Dptrs> vecvecEC;
      for (Int_t iHist = 0; iHist < 4; iHist++) {
         hEC_pos[iHist] =
            new TH1D(ECHistNames[iHist] + "_poschar_" + MpdPidUtils::cParticleShortName[MpdPidUtils::kElectron], "", 50,
                     PMIN, PMAX);
         hEC_pos[iHist]->Sumw2();
         vecEC_pos.push_back(hEC_pos[iHist]);

         hEC_neg[iHist] =
            new TH1D(ECHistNames[iHist] + "_negchar_" + MpdPidUtils::cParticleShortName[MpdPidUtils::kElectron], "", 50,
                     PMIN, PMAX);
         hEC_neg[iHist]->Sumw2();
         vecEC_neg.push_back(hEC_neg[iHist]);
      }
      vecvecEC.push_back(vecEC_pos);
      vecvecEC.push_back(vecEC_neg);
      fEffContMap.insert(pair<MpdPidUtils::ePartType, std::vector<vecTH1Dptrs>>(MpdPidUtils::kElectron, vecvecEC));
   }

   if (Particles.Contains("mu")) {
      /// dE/dx
      if ((fTrackingState == MpdPidUtils::kCF) || (fTrackingState == MpdPidUtils::kCFHM)) {
         MaxEnLossLoc                = 15000.;
         XMIN[0][MpdPidUtils::kMuon] = 2000.;
         XMIN[1][MpdPidUtils::kMuon] = 1000.;
         XMAX[0][MpdPidUtils::kMuon] = 15000.;
         XMAX[1][MpdPidUtils::kMuon] = 8000.;
         for (Int_t i = 2; i < MpdPidUtils::nQAHists; i++) {
            XMIN[i][MpdPidUtils::kMuon] = 1000.;
            XMAX[i][MpdPidUtils::kMuon] = 5000.;
         }
      } else {
         MaxEnLossLoc                = 15.e-06;
         XMAX[0][MpdPidUtils::kMuon] = 15.e-06;
         XMAX[1][MpdPidUtils::kMuon] = 5.e-06;
         for (Int_t i = 0; i < MpdPidUtils::nQAHists; i++) XMIN[i][MpdPidUtils::kMuon] = 1.e-06;
         for (Int_t i = 2; i < MpdPidUtils::nQAHists; i++) XMAX[i][MpdPidUtils::kMuon] = 3.e-06;
      }
      MpdPidUtils::dEdXStruct_t SMu;
      for (Int_t iHist = 0; iHist < MpdPidUtils::nQAHists; iHist++) {
         histName = "Mu_" + TString::Itoa(iHist, 10);
         SMu.dEdXPart[iHist] =
            new TH1D(histName, "", nbins, XMIN[iHist][MpdPidUtils::kMuon], XMAX[iHist][MpdPidUtils::kMuon]);
      }
      SMu.BetheBlochHist = new TH2D("hMuBB_QA", "", 300, 0.05, 3.0, 300, MinEnLoss, MaxEnLossLoc);
      SMu.ibeg           = 0;
      SMu.iend           = 39;
      fEnLossMap.insert(pair<MpdPidUtils::ePartType, MpdPidUtils::dEdXStruct_t>(MpdPidUtils::kMuon, SMu));

      /// m-squared
      TH2D *mass2Mu = new TH2D("mass2Mu", "", nbins, 0., 4., nbins, -0.8, 0.8);
      fMSquaredMap.insert(pair<MpdPidUtils::ePartType, TH2D *>(MpdPidUtils::kMuon, mass2Mu));

      /// Abundance
      vecTH1Dptrs vecYield;
      TH1D       *hYield_mupos = new TH1D("hYield_mupos", "", 600, 0., 3.);
      vecYield.push_back(hYield_mupos);
      TH1D *hYield_muneg = new TH1D("hYield_muneg", "", 600, 0., 3.);
      vecYield.push_back(hYield_muneg);
      fAbundanceMap.insert(pair<MpdPidUtils::ePartType, vecTH1Dptrs>(MpdPidUtils::kMuon, vecYield));

      /// Efficiency & Contamination
      TH1D                    *hEC_pos[4], *hEC_neg[4];
      vecTH1Dptrs              vecEC_pos, vecEC_neg;
      std::vector<vecTH1Dptrs> vecvecEC;
      for (Int_t iHist = 0; iHist < 4; iHist++) {
         hEC_pos[iHist] = new TH1D(
            ECHistNames[iHist] + "_poschar_" + MpdPidUtils::cParticleShortName[MpdPidUtils::kMuon], "", 50, PMIN, PMAX);
         hEC_pos[iHist]->Sumw2();
         vecEC_pos.push_back(hEC_pos[iHist]);

         hEC_neg[iHist] = new TH1D(
            ECHistNames[iHist] + "_negchar_" + MpdPidUtils::cParticleShortName[MpdPidUtils::kMuon], "", 50, PMIN, PMAX);
         hEC_neg[iHist]->Sumw2();
         vecEC_neg.push_back(hEC_neg[iHist]);
      }
      vecvecEC.push_back(vecEC_pos);
      vecvecEC.push_back(vecEC_neg);
      fEffContMap.insert(pair<MpdPidUtils::ePartType, std::vector<vecTH1Dptrs>>(MpdPidUtils::kMuon, vecvecEC));
   }

   if (Particles.Contains("pi")) {
      /// dE/dx
      if (fTrackingState == MpdPidUtils::kCFHM) {
         MaxEnLossLoc = 14.;
         for (Int_t i = 0; i < 9; i++) XMIN[i][MpdPidUtils::kPion] = 0.7;
         for (Int_t i = 9; i < 17; i++) XMIN[i][MpdPidUtils::kPion] = 0.8;
         for (Int_t i = 17; i < 30; i++) XMIN[i][MpdPidUtils::kPion] = 0.9;
         for (Int_t i = 30; i < MpdPidUtils::nQAHists; i++) XMIN[i][MpdPidUtils::kPion] = 1.0;
         XMAX[0][MpdPidUtils::kPion]  = 14.;
         XMAX[1][MpdPidUtils::kPion]  = 6.0;
         XMAX[2][MpdPidUtils::kPion]  = 4.0;
         XMAX[3][MpdPidUtils::kPion]  = 4.0;
         XMAX[4][MpdPidUtils::kPion]  = 3.0;
         XMAX[5][MpdPidUtils::kPion]  = 2.5;
         XMAX[6][MpdPidUtils::kPion]  = 2.5;
         XMAX[16][MpdPidUtils::kPion] = 2.2;
         for (Int_t i = 7; i < 13; i++) XMAX[i][MpdPidUtils::kPion] = 2.4;
         for (Int_t i = 13; i < 16; i++) XMAX[i][MpdPidUtils::kPion] = 2.3;
         for (Int_t i = 17; i < 35; i++) XMAX[i][MpdPidUtils::kPion] = 2.1;
         for (Int_t i = 35; i < 38; i++) XMAX[i][MpdPidUtils::kPion] = 2.2;
         for (Int_t i = 38; i < MpdPidUtils::nQAHists; i++) XMAX[i][MpdPidUtils::kPion] = 2.3;
      } else if (fTrackingState == MpdPidUtils::kCF) {
         MaxEnLossLoc                = 25000.;
         XMIN[0][MpdPidUtils::kPion] = 1000.;
         XMAX[0][MpdPidUtils::kPion] = 25000.;
         XMAX[1][MpdPidUtils::kPion] = 10000.;
         XMAX[2][MpdPidUtils::kPion] = 8000.;
         for (Int_t i = 1; i < MpdPidUtils::nQAHists; i++) XMIN[i][MpdPidUtils::kPion] = 2000.;
         for (Int_t i = 3; i < MpdPidUtils::nQAHists; i++) XMAX[i][MpdPidUtils::kPion] = 6000.;
      } else {
         MaxEnLossLoc                = 15.e-06;
         XMAX[0][MpdPidUtils::kPion] = 15.e-06;
         XMAX[1][MpdPidUtils::kPion] = 6.e-06;
         XMAX[2][MpdPidUtils::kPion] = 3.5e-06;
         XMAX[3][MpdPidUtils::kPion] = 3.e-06;
         XMAX[4][MpdPidUtils::kPion] = 2.5e-06;
         XMAX[5][MpdPidUtils::kPion] = 2.5e-06;
         XMAX[6][MpdPidUtils::kPion] = 2.2e-06;
         for (Int_t i = 0; i < MpdPidUtils::nQAHists; i++) XMIN[i][MpdPidUtils::kPion] = 1.e-06;
         for (Int_t i = 7; i < MpdPidUtils::nQAHists; i++) XMAX[i][MpdPidUtils::kPion] = 2.1e-06;
      }
      MpdPidUtils::dEdXStruct_t SPi;
      for (Int_t iHist = 0; iHist < MpdPidUtils::nQAHists; iHist++) {
         histName = "Pi_" + TString::Itoa(iHist, 10);
         SPi.dEdXPart[iHist] =
            new TH1D(histName, "", nbins, XMIN[iHist][MpdPidUtils::kPion], XMAX[iHist][MpdPidUtils::kPion]);
      }
      SPi.BetheBlochHist = new TH2D("hPiBB_QA", "", 300, 0.05, 3.0, 300, MinEnLoss, MaxEnLossLoc);
      SPi.ibeg           = 0;
      SPi.iend           = 39;
      fEnLossMap.insert(pair<MpdPidUtils::ePartType, MpdPidUtils::dEdXStruct_t>(MpdPidUtils::kPion, SPi));

      /// m-squared
      TH2D *mass2Pi = new TH2D("mass2Pi", "", nbins, 0., 4., nbins, -0.5, 1.5);
      fMSquaredMap.insert(pair<MpdPidUtils::ePartType, TH2D *>(MpdPidUtils::kPion, mass2Pi));

      /// Abundance
      vecTH1Dptrs vecYield;
      TH1D       *hYield_pipos = new TH1D("hYield_pipos", "", 600, 0., 3.);
      vecYield.push_back(hYield_pipos);
      TH1D *hYield_pineg = new TH1D("hYield_pineg", "", 600, 0., 3.);
      vecYield.push_back(hYield_pineg);
      fAbundanceMap.insert(pair<MpdPidUtils::ePartType, vecTH1Dptrs>(MpdPidUtils::kPion, vecYield));

      /// Efficiency & Contamination
      TH1D                    *hEC_pos[4], *hEC_neg[4];
      vecTH1Dptrs              vecEC_pos, vecEC_neg;
      std::vector<vecTH1Dptrs> vecvecEC;
      for (Int_t iHist = 0; iHist < 4; iHist++) {
         hEC_pos[iHist] = new TH1D(
            ECHistNames[iHist] + "_poschar_" + MpdPidUtils::cParticleShortName[MpdPidUtils::kPion], "", 50, PMIN, PMAX);
         hEC_pos[iHist]->Sumw2();
         vecEC_pos.push_back(hEC_pos[iHist]);

         hEC_neg[iHist] = new TH1D(
            ECHistNames[iHist] + "_negchar_" + MpdPidUtils::cParticleShortName[MpdPidUtils::kPion], "", 50, PMIN, PMAX);
         hEC_neg[iHist]->Sumw2();
         vecEC_neg.push_back(hEC_neg[iHist]);
      }
      vecvecEC.push_back(vecEC_pos);
      vecvecEC.push_back(vecEC_neg);
      fEffContMap.insert(pair<MpdPidUtils::ePartType, std::vector<vecTH1Dptrs>>(MpdPidUtils::kPion, vecvecEC));
   }

   if (Particles.Contains("ka")) {
      /// dE/dx
      if (fTrackingState == MpdPidUtils::kCFHM) {
         MaxEnLossLoc                 = 70;
         XMIN[0][MpdPidUtils::kKaon]  = 10.;
         XMIN[1][MpdPidUtils::kKaon]  = 6.;
         XMIN[2][MpdPidUtils::kKaon]  = 5.;
         XMIN[3][MpdPidUtils::kKaon]  = 3.;
         XMIN[4][MpdPidUtils::kKaon]  = 2.;
         XMIN[5][MpdPidUtils::kKaon]  = 2.;
         XMIN[21][MpdPidUtils::kKaon] = 0.9;
         XMIN[22][MpdPidUtils::kKaon] = 0.9;
         for (Int_t i = 6; i < 21; i++) XMIN[i][MpdPidUtils::kKaon] = 1.;
         for (Int_t i = 23; i < MpdPidUtils::nQAHists; i++) XMIN[i][MpdPidUtils::kKaon] = 0.8;
         XMAX[0][MpdPidUtils::kKaon]  = 70.;
         XMAX[1][MpdPidUtils::kKaon]  = 40.;
         XMAX[2][MpdPidUtils::kKaon]  = 20.;
         XMAX[3][MpdPidUtils::kKaon]  = 12.;
         XMAX[4][MpdPidUtils::kKaon]  = 8.;
         XMAX[5][MpdPidUtils::kKaon]  = 6.;
         XMAX[6][MpdPidUtils::kKaon]  = 5.5;
         XMAX[7][MpdPidUtils::kKaon]  = 4.5;
         XMAX[8][MpdPidUtils::kKaon]  = 4.;
         XMAX[9][MpdPidUtils::kKaon]  = 3.5;
         XMAX[10][MpdPidUtils::kKaon] = 3.5;
         XMAX[11][MpdPidUtils::kKaon] = 3.;
         XMAX[12][MpdPidUtils::kKaon] = 2.8;
         XMAX[13][MpdPidUtils::kKaon] = 2.6;
         XMAX[14][MpdPidUtils::kKaon] = 2.4;
         XMAX[15][MpdPidUtils::kKaon] = 2.3;
         XMAX[16][MpdPidUtils::kKaon] = 2.2;
         XMAX[17][MpdPidUtils::kKaon] = 2.2;
         XMAX[18][MpdPidUtils::kKaon] = 2.1;
         for (Int_t i = 19; i < 22; i++) XMAX[i][MpdPidUtils::kKaon] = 2.;
         for (Int_t i = 22; i < 25; i++) XMAX[i][MpdPidUtils::kKaon] = 1.9;
         for (Int_t i = 25; i < 38; i++) XMAX[i][MpdPidUtils::kKaon] = 1.8;
         for (Int_t i = 38; i < MpdPidUtils::nQAHists; i++) XMAX[i][MpdPidUtils::kKaon] = 1.7;
      } else if (fTrackingState == MpdPidUtils::kCF) {
         MaxEnLossLoc                = 200000.;
         XMIN[0][MpdPidUtils::kKaon] = 10000.;
         XMIN[1][MpdPidUtils::kKaon] = 12000.;
         XMIN[2][MpdPidUtils::kKaon] = 8000.;
         XMIN[3][MpdPidUtils::kKaon] = 5000.;
         XMIN[4][MpdPidUtils::kKaon] = 3000.;
         for (Int_t i = 5; i < 14; i++) XMIN[i][MpdPidUtils::kKaon] = 2000.;
         for (Int_t i = 14; i < MpdPidUtils::nQAHists; i++) XMIN[i][MpdPidUtils::kKaon] = 1000.;
         XMAX[0][MpdPidUtils::kKaon]  = 100000.;
         XMAX[1][MpdPidUtils::kKaon]  = 60000.;
         XMAX[2][MpdPidUtils::kKaon]  = 40000.;
         XMAX[3][MpdPidUtils::kKaon]  = 22000.;
         XMAX[4][MpdPidUtils::kKaon]  = 17000.;
         XMAX[5][MpdPidUtils::kKaon]  = 13000.;
         XMAX[6][MpdPidUtils::kKaon]  = 13000.;
         XMAX[7][MpdPidUtils::kKaon]  = 10000.;
         XMAX[8][MpdPidUtils::kKaon]  = 10000.;
         XMAX[9][MpdPidUtils::kKaon]  = 9000.;
         XMAX[10][MpdPidUtils::kKaon] = 9000.;
         XMAX[11][MpdPidUtils::kKaon] = 7000.;
         XMAX[12][MpdPidUtils::kKaon] = 7000.;
         XMAX[13][MpdPidUtils::kKaon] = 6000.;
         XMAX[14][MpdPidUtils::kKaon] = 6000.;
         for (Int_t i = 15; i < MpdPidUtils::nQAHists; i++) XMAX[i][MpdPidUtils::kKaon] = 5000.;
      } else {
         MaxEnLossLoc                 = 100.e-06;
         XMIN[0][MpdPidUtils::kKaon]  = 10.e-06;
         XMIN[1][MpdPidUtils::kKaon]  = 10.e-06;
         XMIN[2][MpdPidUtils::kKaon]  = 3.e-06;
         XMIN[3][MpdPidUtils::kKaon]  = 3.e-06;
         XMIN[4][MpdPidUtils::kKaon]  = 3.e-06;
         XMIN[5][MpdPidUtils::kKaon]  = 2.e-06;
         XMIN[6][MpdPidUtils::kKaon]  = 2.e-06;
         XMIN[7][MpdPidUtils::kKaon]  = 1.5e-06;
         XMIN[8][MpdPidUtils::kKaon]  = 1.5e-06;
         XMIN[9][MpdPidUtils::kKaon]  = 1.5e-06;
         XMIN[10][MpdPidUtils::kKaon] = 1.5e-06;
         for (Int_t i = 11; i < MpdPidUtils::nQAHists; i++) XMIN[i][MpdPidUtils::kKaon] = 1.e-06;
         XMAX[0][MpdPidUtils::kKaon]  = 100.e-06;
         XMAX[1][MpdPidUtils::kKaon]  = 50.e-06;
         XMAX[2][MpdPidUtils::kKaon]  = 25.e-06;
         XMAX[3][MpdPidUtils::kKaon]  = 12.e-06;
         XMAX[4][MpdPidUtils::kKaon]  = 8.e-06;
         XMAX[5][MpdPidUtils::kKaon]  = 7.e-06;
         XMAX[6][MpdPidUtils::kKaon]  = 5.e-06;
         XMAX[7][MpdPidUtils::kKaon]  = 4.5e-06;
         XMAX[8][MpdPidUtils::kKaon]  = 4.e-06;
         XMAX[9][MpdPidUtils::kKaon]  = 3.5e-06;
         XMAX[10][MpdPidUtils::kKaon] = 3.e-06;
         XMAX[11][MpdPidUtils::kKaon] = 3.e-06;
         XMAX[12][MpdPidUtils::kKaon] = 3.e-06;
         XMAX[13][MpdPidUtils::kKaon] = 3.e-06;
         XMAX[14][MpdPidUtils::kKaon] = 3.e-06;
         for (Int_t i = 15; i < MpdPidUtils::nQAHists; i++) XMAX[i][MpdPidUtils::kKaon] = 2.5e-06;
      }
      MpdPidUtils::dEdXStruct_t SKa;
      for (Int_t iHist = 0; iHist < MpdPidUtils::nQAHists; iHist++) {
         histName = "Ka_" + TString::Itoa(iHist, 10);
         SKa.dEdXPart[iHist] =
            new TH1D(histName, "", nbins, XMIN[iHist][MpdPidUtils::kKaon], XMAX[iHist][MpdPidUtils::kKaon]);
      }
      SKa.BetheBlochHist = new TH2D("hKaBB_QA", "", 300, 0.05, 3.0, 300, MinEnLoss, MaxEnLossLoc);
      SKa.ibeg           = 0;
      SKa.iend           = 39;
      fEnLossMap.insert(pair<MpdPidUtils::ePartType, MpdPidUtils::dEdXStruct_t>(MpdPidUtils::kKaon, SKa));

      /// m-squared
      TH2D *mass2Ka = new TH2D("mass2Ka", "", nbins, 0., 4., nbins, -0.5, 1.5);
      fMSquaredMap.insert(pair<MpdPidUtils::ePartType, TH2D *>(MpdPidUtils::kKaon, mass2Ka));

      /// Abundance
      vecTH1Dptrs vecYield;
      TH1D       *hYield_kapos = new TH1D("hYield_kapos", "", 600, 0., 3.);
      vecYield.push_back(hYield_kapos);
      TH1D *hYield_kaneg = new TH1D("hYield_kaneg", "", 600, 0., 3.);
      vecYield.push_back(hYield_kaneg);
      fAbundanceMap.insert(pair<MpdPidUtils::ePartType, vecTH1Dptrs>(MpdPidUtils::kKaon, vecYield));

      /// Efficiency & Contamination
      TH1D                    *hEC_pos[4], *hEC_neg[4];
      vecTH1Dptrs              vecEC_pos, vecEC_neg;
      std::vector<vecTH1Dptrs> vecvecEC;
      for (Int_t iHist = 0; iHist < 4; iHist++) {
         hEC_pos[iHist] = new TH1D(
            ECHistNames[iHist] + "_poschar_" + MpdPidUtils::cParticleShortName[MpdPidUtils::kKaon], "", 50, PMIN, PMAX);
         hEC_pos[iHist]->Sumw2();
         vecEC_pos.push_back(hEC_pos[iHist]);

         hEC_neg[iHist] = new TH1D(
            ECHistNames[iHist] + "_negchar_" + MpdPidUtils::cParticleShortName[MpdPidUtils::kKaon], "", 50, PMIN, PMAX);
         hEC_neg[iHist]->Sumw2();
         vecEC_neg.push_back(hEC_neg[iHist]);
      }
      vecvecEC.push_back(vecEC_pos);
      vecvecEC.push_back(vecEC_neg);
      fEffContMap.insert(pair<MpdPidUtils::ePartType, std::vector<vecTH1Dptrs>>(MpdPidUtils::kKaon, vecvecEC));
   }

   if (Particles.Contains("pr")) {
      /// dE/dx
      if (fTrackingState == MpdPidUtils::kCFHM) {
         MaxEnLossLoc                   = 90.;
         XMIN[0][MpdPidUtils::kProton]  = 30.;
         XMIN[1][MpdPidUtils::kProton]  = 12.;
         XMIN[2][MpdPidUtils::kProton]  = 10.;
         XMIN[3][MpdPidUtils::kProton]  = 10.;
         XMIN[4][MpdPidUtils::kProton]  = 6.;
         XMIN[5][MpdPidUtils::kProton]  = 6.;
         XMIN[6][MpdPidUtils::kProton]  = 5.;
         XMIN[7][MpdPidUtils::kProton]  = 4.;
         XMIN[8][MpdPidUtils::kProton]  = 3.5;
         XMIN[9][MpdPidUtils::kProton]  = 3.;
         XMIN[10][MpdPidUtils::kProton] = 2.7;
         XMIN[11][MpdPidUtils::kProton] = 2.5;
         XMIN[12][MpdPidUtils::kProton] = 2.2;
         XMIN[13][MpdPidUtils::kProton] = 2.;
         XMIN[14][MpdPidUtils::kProton] = 1.8;
         XMIN[15][MpdPidUtils::kProton] = 1.7;
         XMIN[16][MpdPidUtils::kProton] = 1.6;
         XMIN[17][MpdPidUtils::kProton] = 1.5;
         XMIN[18][MpdPidUtils::kProton] = 1.4;
         XMIN[19][MpdPidUtils::kProton] = 1.4;
         for (Int_t i = 20; i < 23; i++) XMIN[i][MpdPidUtils::kProton] = 1.3;
         for (Int_t i = 23; i < 26; i++) XMIN[i][MpdPidUtils::kProton] = 1.2;
         for (Int_t i = 26; i < 31; i++) XMIN[i][MpdPidUtils::kProton] = 1.1;
         for (Int_t i = 31; i < 35; i++) XMIN[i][MpdPidUtils::kProton] = 1.0;
         for (Int_t i = 35; i < MpdPidUtils::nQAHists; i++) XMIN[i][MpdPidUtils::kProton] = 0.9;
         XMAX[0][MpdPidUtils::kProton]  = 90.;
         XMAX[1][MpdPidUtils::kProton]  = 80.;
         XMAX[2][MpdPidUtils::kProton]  = 60.;
         XMAX[3][MpdPidUtils::kProton]  = 37.;
         XMAX[4][MpdPidUtils::kProton]  = 25.;
         XMAX[5][MpdPidUtils::kProton]  = 16.;
         XMAX[6][MpdPidUtils::kProton]  = 12.;
         XMAX[7][MpdPidUtils::kProton]  = 10.;
         XMAX[8][MpdPidUtils::kProton]  = 8.;
         XMAX[9][MpdPidUtils::kProton]  = 7.;
         XMAX[10][MpdPidUtils::kProton] = 5.8;
         XMAX[11][MpdPidUtils::kProton] = 5.2;
         XMAX[12][MpdPidUtils::kProton] = 4.6;
         XMAX[13][MpdPidUtils::kProton] = 4.2;
         XMAX[14][MpdPidUtils::kProton] = 3.8;
         XMAX[15][MpdPidUtils::kProton] = 3.5;
         XMAX[16][MpdPidUtils::kProton] = 3.4;
         XMAX[17][MpdPidUtils::kProton] = 3.3;
         XMAX[18][MpdPidUtils::kProton] = 3.2;
         XMAX[19][MpdPidUtils::kProton] = 3.0;
         XMAX[20][MpdPidUtils::kProton] = 2.8;
         XMAX[21][MpdPidUtils::kProton] = 2.7;
         XMAX[22][MpdPidUtils::kProton] = 2.6;
         XMAX[23][MpdPidUtils::kProton] = 2.5;
         XMAX[24][MpdPidUtils::kProton] = 2.5;
         XMAX[25][MpdPidUtils::kProton] = 2.5;
         XMAX[26][MpdPidUtils::kProton] = 2.4;
         XMAX[27][MpdPidUtils::kProton] = 2.3;
         XMAX[28][MpdPidUtils::kProton] = 2.3;
         XMAX[29][MpdPidUtils::kProton] = 2.2;
         XMAX[30][MpdPidUtils::kProton] = 2.2;
         XMAX[31][MpdPidUtils::kProton] = 2.1;
         XMAX[32][MpdPidUtils::kProton] = 2.1;
         for (Int_t i = 33; i < MpdPidUtils::nQAHists; i++) XMAX[i][MpdPidUtils::kProton] = 2.0;
      } else if (fTrackingState == MpdPidUtils::kCF) {
         MaxEnLossLoc                  = 400000.;
         XMIN[0][MpdPidUtils::kProton] = 100000.;
         XMIN[1][MpdPidUtils::kProton] = 50000.;
         XMIN[2][MpdPidUtils::kProton] = 20000.;
         XMIN[3][MpdPidUtils::kProton] = 15000.;
         XMIN[4][MpdPidUtils::kProton] = 10000.;
         XMIN[5][MpdPidUtils::kProton] = 10000.;
         XMIN[6][MpdPidUtils::kProton] = 8000.;
         XMIN[7][MpdPidUtils::kProton] = 8000.;
         XMIN[8][MpdPidUtils::kProton] = 8000.;
         for (Int_t i = 9; i < 14; i++) XMIN[i][MpdPidUtils::kProton] = 5000.;
         XMIN[14][MpdPidUtils::kProton] = 2000.;
         for (Int_t i = 15; i < 23; i++) XMIN[i][MpdPidUtils::kProton] = 3000.;
         for (Int_t i = 23; i < MpdPidUtils::nQAHists; i++) XMIN[i][MpdPidUtils::kProton] = 2000.;
         XMAX[0][MpdPidUtils::kProton]  = 800000.;
         XMAX[1][MpdPidUtils::kProton]  = 250000.;
         XMAX[2][MpdPidUtils::kProton]  = 180000.;
         XMAX[3][MpdPidUtils::kProton]  = 90000.;
         XMAX[4][MpdPidUtils::kProton]  = 50000.;
         XMAX[5][MpdPidUtils::kProton]  = 30000.;
         XMAX[6][MpdPidUtils::kProton]  = 25000.;
         XMAX[7][MpdPidUtils::kProton]  = 19000.;
         XMAX[8][MpdPidUtils::kProton]  = 17000.;
         XMAX[9][MpdPidUtils::kProton]  = 15000.;
         XMAX[10][MpdPidUtils::kProton] = 15000.;
         XMAX[11][MpdPidUtils::kProton] = 15000.;
         XMAX[12][MpdPidUtils::kProton] = 12000.;
         XMAX[13][MpdPidUtils::kProton] = 12000.;
         XMAX[14][MpdPidUtils::kProton] = 12000.;
         XMAX[15][MpdPidUtils::kProton] = 9000.;
         XMAX[16][MpdPidUtils::kProton] = 9000.;
         XMAX[17][MpdPidUtils::kProton] = 8000.;
         XMAX[18][MpdPidUtils::kProton] = 8000.;
         XMAX[19][MpdPidUtils::kProton] = 8000.;
         XMAX[20][MpdPidUtils::kProton] = 7000.;
         XMAX[21][MpdPidUtils::kProton] = 7000.;
         XMAX[22][MpdPidUtils::kProton] = 7000.;
         for (Int_t i = 23; i < MpdPidUtils::nQAHists; i++) XMAX[i][MpdPidUtils::kProton] = 6000.;
      } else {
         MaxEnLossLoc                   = 180.e-06;
         XMIN[0][MpdPidUtils::kProton]  = 75.e-06;
         XMIN[1][MpdPidUtils::kProton]  = 40.e-06;
         XMIN[2][MpdPidUtils::kProton]  = 20.e-06;
         XMIN[3][MpdPidUtils::kProton]  = 10.e-06;
         XMIN[4][MpdPidUtils::kProton]  = 10.e-06;
         XMIN[5][MpdPidUtils::kProton]  = 8.e-06;
         XMIN[6][MpdPidUtils::kProton]  = 6.e-06;
         XMIN[7][MpdPidUtils::kProton]  = 5.e-06;
         XMIN[8][MpdPidUtils::kProton]  = 4.e-06;
         XMIN[9][MpdPidUtils::kProton]  = 3.5e-06;
         XMIN[10][MpdPidUtils::kProton] = 3.e-06;
         XMIN[11][MpdPidUtils::kProton] = 3.e-06;
         XMIN[12][MpdPidUtils::kProton] = 2.5e-06;
         for (Int_t i = 13; i < 17; i++) XMIN[i][MpdPidUtils::kProton] = 2.e-06;
         for (Int_t i = 17; i < 23; i++) XMIN[i][MpdPidUtils::kProton] = 1.5e-06;
         for (Int_t i = 23; i < MpdPidUtils::nQAHists; i++) XMIN[i][MpdPidUtils::kProton] = 1.e-06;
         XMAX[0][MpdPidUtils::kProton]  = 180.e-06;
         XMAX[1][MpdPidUtils::kProton]  = 100.e-06;
         XMAX[2][MpdPidUtils::kProton]  = 70.e-06;
         XMAX[3][MpdPidUtils::kProton]  = 35.e-06;
         XMAX[4][MpdPidUtils::kProton]  = 25.e-06;
         XMAX[5][MpdPidUtils::kProton]  = 17.e-06;
         XMAX[6][MpdPidUtils::kProton]  = 13.e-06;
         XMAX[7][MpdPidUtils::kProton]  = 11.e-06;
         XMAX[8][MpdPidUtils::kProton]  = 9.e-06;
         XMAX[9][MpdPidUtils::kProton]  = 7.5e-06;
         XMAX[10][MpdPidUtils::kProton] = 6.5e-06;
         XMAX[11][MpdPidUtils::kProton] = 6.e-06;
         XMAX[12][MpdPidUtils::kProton] = 5.e-06;
         XMAX[13][MpdPidUtils::kProton] = 5.e-06;
         XMAX[14][MpdPidUtils::kProton] = 4.5e-06;
         XMAX[15][MpdPidUtils::kProton] = 4.e-06;
         XMAX[16][MpdPidUtils::kProton] = 4.e-06;
         for (Int_t i = 17; i < 21; i++) XMAX[i][MpdPidUtils::kProton] = 3.5e-06;
         for (Int_t i = 21; i < 31; i++) XMAX[i][MpdPidUtils::kProton] = 3.e-06;
         for (Int_t i = 31; i < MpdPidUtils::nQAHists; i++) XMAX[i][MpdPidUtils::kProton] = 2.5e-06;
      }
      MpdPidUtils::dEdXStruct_t SPr;
      for (Int_t iHist = 0; iHist < MpdPidUtils::nQAHists; iHist++) {
         histName = "Pr_" + TString::Itoa(iHist, 10);
         SPr.dEdXPart[iHist] =
            new TH1D(histName, "", nbins, XMIN[iHist][MpdPidUtils::kProton], XMAX[iHist][MpdPidUtils::kProton]);
      }
      SPr.BetheBlochHist = new TH2D("hPrBB_QA", "", 300, 0.05, 3.0, 300, MinEnLoss, MaxEnLossLoc);
      SPr.ibeg           = 0;
      SPr.iend           = 39;
      fEnLossMap.insert(pair<MpdPidUtils::ePartType, MpdPidUtils::dEdXStruct_t>(MpdPidUtils::kProton, SPr));

      /// m-squared
      TH2D *mass2Pr = new TH2D("mass2Pr", "", nbins, 0., 4., nbins, -0.5, 1.5);
      fMSquaredMap.insert(pair<MpdPidUtils::ePartType, TH2D *>(MpdPidUtils::kProton, mass2Pr));

      /// Abundance
      vecTH1Dptrs vecYield;
      TH1D       *hYield_prpos = new TH1D("hYield_prpos", "", 600, 0., 3.);
      vecYield.push_back(hYield_prpos);
      TH1D *hYield_prneg = new TH1D("hYield_prneg", "", 600, 0., 3.);
      vecYield.push_back(hYield_prneg);
      fAbundanceMap.insert(pair<MpdPidUtils::ePartType, vecTH1Dptrs>(MpdPidUtils::kProton, vecYield));

      /// Efficiency & Contamination
      TH1D                    *hEC_pos[4], *hEC_neg[4];
      vecTH1Dptrs              vecEC_pos, vecEC_neg;
      std::vector<vecTH1Dptrs> vecvecEC;
      for (Int_t iHist = 0; iHist < 4; iHist++) {
         hEC_pos[iHist] =
            new TH1D(ECHistNames[iHist] + "_poschar_" + MpdPidUtils::cParticleShortName[MpdPidUtils::kProton], "", 50,
                     PMIN, PMAX);
         hEC_pos[iHist]->Sumw2();
         vecEC_pos.push_back(hEC_pos[iHist]);

         hEC_neg[iHist] =
            new TH1D(ECHistNames[iHist] + "_negchar_" + MpdPidUtils::cParticleShortName[MpdPidUtils::kProton], "", 50,
                     PMIN, PMAX);
         hEC_neg[iHist]->Sumw2();
         vecEC_neg.push_back(hEC_neg[iHist]);
      }
      vecvecEC.push_back(vecEC_pos);
      vecvecEC.push_back(vecEC_neg);
      fEffContMap.insert(pair<MpdPidUtils::ePartType, std::vector<vecTH1Dptrs>>(MpdPidUtils::kProton, vecvecEC));
   }

   if (Particles.Contains("de")) {
      /// dE/dx
      if (fTrackingState == MpdPidUtils::kCFHM) {
         MaxEnLossLoc                     = 90.;
         XMIN[0][MpdPidUtils::kDeuteron]  = 20.;
         XMIN[1][MpdPidUtils::kDeuteron]  = 20.;
         XMIN[2][MpdPidUtils::kDeuteron]  = 15.;
         XMIN[3][MpdPidUtils::kDeuteron]  = 15.;
         XMIN[4][MpdPidUtils::kDeuteron]  = 10.;
         XMIN[5][MpdPidUtils::kDeuteron]  = 10.;
         XMIN[6][MpdPidUtils::kDeuteron]  = 8.;
         XMIN[7][MpdPidUtils::kDeuteron]  = 7.;
         XMIN[8][MpdPidUtils::kDeuteron]  = 6.;
         XMIN[9][MpdPidUtils::kDeuteron]  = 5.;
         XMIN[10][MpdPidUtils::kDeuteron] = 5.;
         XMIN[11][MpdPidUtils::kDeuteron] = 5.;
         XMIN[12][MpdPidUtils::kDeuteron] = 4.;
         XMIN[13][MpdPidUtils::kDeuteron] = 3.;
         XMIN[14][MpdPidUtils::kDeuteron] = 3.;
         XMIN[15][MpdPidUtils::kDeuteron] = 2.;
         XMIN[16][MpdPidUtils::kDeuteron] = 2.;
         XMIN[17][MpdPidUtils::kDeuteron] = 1.5;
         XMIN[18][MpdPidUtils::kDeuteron] = 1.5;
         XMIN[19][MpdPidUtils::kDeuteron] = 1.5;
         for (Int_t i = 20; i < 30; i++) XMIN[i][MpdPidUtils::kDeuteron] = 1.0;
         for (Int_t i = 30; i < 37; i++) XMIN[i][MpdPidUtils::kDeuteron] = 0.5;
         XMAX[0][MpdPidUtils::kDeuteron]  = 90.;
         XMAX[1][MpdPidUtils::kDeuteron]  = 70.;
         XMAX[2][MpdPidUtils::kDeuteron]  = 55.;
         XMAX[3][MpdPidUtils::kDeuteron]  = 45.;
         XMAX[4][MpdPidUtils::kDeuteron]  = 40.;
         XMAX[5][MpdPidUtils::kDeuteron]  = 32.;
         XMAX[6][MpdPidUtils::kDeuteron]  = 26.;
         XMAX[7][MpdPidUtils::kDeuteron]  = 22.;
         XMAX[8][MpdPidUtils::kDeuteron]  = 19.;
         XMAX[9][MpdPidUtils::kDeuteron]  = 16.;
         XMAX[10][MpdPidUtils::kDeuteron] = 14.;
         XMAX[11][MpdPidUtils::kDeuteron] = 12.;
         XMAX[12][MpdPidUtils::kDeuteron] = 11.;
         XMAX[13][MpdPidUtils::kDeuteron] = 10.;
         XMAX[14][MpdPidUtils::kDeuteron] = 9.;
         XMAX[15][MpdPidUtils::kDeuteron] = 9.;
         XMAX[16][MpdPidUtils::kDeuteron] = 8.;
         XMAX[17][MpdPidUtils::kDeuteron] = 7.5;
         XMAX[18][MpdPidUtils::kDeuteron] = 7.5;
         XMAX[19][MpdPidUtils::kDeuteron] = 7.;
         XMAX[20][MpdPidUtils::kDeuteron] = 6.5;
         XMAX[21][MpdPidUtils::kDeuteron] = 6.5;
         XMAX[22][MpdPidUtils::kDeuteron] = 6.;
         XMAX[23][MpdPidUtils::kDeuteron] = 5.5;
         XMAX[24][MpdPidUtils::kDeuteron] = 5.5;
         XMAX[25][MpdPidUtils::kDeuteron] = 5.;
         XMAX[26][MpdPidUtils::kDeuteron] = 5.;
         XMAX[27][MpdPidUtils::kDeuteron] = 5.;
         XMAX[28][MpdPidUtils::kDeuteron] = 4.5;
         for (Int_t i = 29; i < 35; i++) XMAX[i][MpdPidUtils::kDeuteron] = 4.;
         XMAX[35][MpdPidUtils::kDeuteron] = 3.5;
         XMAX[36][MpdPidUtils::kDeuteron] = 3.5;
      } else if (fTrackingState == MpdPidUtils::kCF) {
         MaxEnLossLoc                    = 500000.;
         XMIN[0][MpdPidUtils::kDeuteron] = 50000.;
         XMIN[1][MpdPidUtils::kDeuteron] = 20000.;
         XMIN[2][MpdPidUtils::kDeuteron] = 10000.;
         XMIN[3][MpdPidUtils::kDeuteron] = 8000.;
         XMIN[4][MpdPidUtils::kDeuteron] = 8000.;
         XMIN[5][MpdPidUtils::kDeuteron] = 5000.;
         for (Int_t i = 6; i < (MpdPidUtils::nQAHists - 2); i++) XMIN[i][MpdPidUtils::kDeuteron] = 2000.;
         XMAX[0][MpdPidUtils::kDeuteron]  = 500000.;
         XMAX[1][MpdPidUtils::kDeuteron]  = 400000.;
         XMAX[2][MpdPidUtils::kDeuteron]  = 250000.;
         XMAX[3][MpdPidUtils::kDeuteron]  = 200000.;
         XMAX[4][MpdPidUtils::kDeuteron]  = 120000.;
         XMAX[5][MpdPidUtils::kDeuteron]  = 80000.;
         XMAX[6][MpdPidUtils::kDeuteron]  = 70000.;
         XMAX[7][MpdPidUtils::kDeuteron]  = 60000.;
         XMAX[8][MpdPidUtils::kDeuteron]  = 50000.;
         XMAX[9][MpdPidUtils::kDeuteron]  = 40000.;
         XMAX[10][MpdPidUtils::kDeuteron] = 35000.;
         XMAX[11][MpdPidUtils::kDeuteron] = 30000.;
         XMAX[12][MpdPidUtils::kDeuteron] = 25000.;
         XMAX[13][MpdPidUtils::kDeuteron] = 20000.;
         XMAX[14][MpdPidUtils::kDeuteron] = 20000.;
         XMAX[15][MpdPidUtils::kDeuteron] = 18000.;
         XMAX[16][MpdPidUtils::kDeuteron] = 15000.;
         XMAX[17][MpdPidUtils::kDeuteron] = 15000.;
         XMAX[18][MpdPidUtils::kDeuteron] = 15000.;
         XMAX[19][MpdPidUtils::kDeuteron] = 13000.;
         XMAX[20][MpdPidUtils::kDeuteron] = 13000.;
         XMAX[21][MpdPidUtils::kDeuteron] = 13000.;
         XMAX[22][MpdPidUtils::kDeuteron] = 13000.;
         XMAX[23][MpdPidUtils::kDeuteron] = 11000.;
         XMAX[24][MpdPidUtils::kDeuteron] = 11000.;
         XMAX[25][MpdPidUtils::kDeuteron] = 10000.;
         XMAX[26][MpdPidUtils::kDeuteron] = 10000.;
         XMAX[27][MpdPidUtils::kDeuteron] = 10000.;
         XMAX[28][MpdPidUtils::kDeuteron] = 9000.;
         XMAX[29][MpdPidUtils::kDeuteron] = 8000.;
         XMAX[30][MpdPidUtils::kDeuteron] = 8000.;
         XMAX[31][MpdPidUtils::kDeuteron] = 7000.;
         XMAX[32][MpdPidUtils::kDeuteron] = 7000.;
         XMAX[33][MpdPidUtils::kDeuteron] = 7000.;
         XMAX[34][MpdPidUtils::kDeuteron] = 7000.;
         XMAX[35][MpdPidUtils::kDeuteron] = 6000.;
         XMAX[36][MpdPidUtils::kDeuteron] = 6000.;
         XMAX[37][MpdPidUtils::kDeuteron] = 6000.;
      } else {
         MaxEnLossLoc = 200.e-06;
         for (Int_t i = 0; i < (MpdPidUtils::nQAHists - 2); i++) XMIN[i][MpdPidUtils::kDeuteron] = 1.e-06;
         XMAX[0][MpdPidUtils::kDeuteron]  = 200.e-06;
         XMAX[1][MpdPidUtils::kDeuteron]  = 200.e-06;
         XMAX[2][MpdPidUtils::kDeuteron]  = 150.e-06;
         XMAX[3][MpdPidUtils::kDeuteron]  = 140.e-06;
         XMAX[4][MpdPidUtils::kDeuteron]  = 100.e-06;
         XMAX[5][MpdPidUtils::kDeuteron]  = 80.e-06;
         XMAX[6][MpdPidUtils::kDeuteron]  = 50.e-06;
         XMAX[7][MpdPidUtils::kDeuteron]  = 40.e-06;
         XMAX[8][MpdPidUtils::kDeuteron]  = 30.e-06;
         XMAX[9][MpdPidUtils::kDeuteron]  = 20.e-06;
         XMAX[10][MpdPidUtils::kDeuteron] = 10.e-06;
         XMAX[11][MpdPidUtils::kDeuteron] = 8.e-06;
         for (Int_t i = 12; i < (MpdPidUtils::nQAHists - 2); i++) XMAX[i][MpdPidUtils::kDeuteron] = 5.e-06;
      }
      MpdPidUtils::dEdXStruct_t SDe;
      for (Int_t iHist = 0; iHist < MpdPidUtils::nQAHists; iHist++) {
         histName = "De_" + TString::Itoa(iHist, 10);
         SDe.dEdXPart[iHist] =
            new TH1D(histName, "", nbins, XMIN[iHist][MpdPidUtils::kDeuteron], XMAX[iHist][MpdPidUtils::kDeuteron]);
      }
      SDe.BetheBlochHist = new TH2D("hDeBB_QA", "", 300, 0.05, 3.0, 300, MinEnLoss, MaxEnLossLoc);
      SDe.ibeg           = (fTrackingState == MpdPidUtils::kCFHM) ? 3 : 2;
      SDe.iend           = 39;
      fEnLossMap.insert(pair<MpdPidUtils::ePartType, MpdPidUtils::dEdXStruct_t>(MpdPidUtils::kDeuteron, SDe));

      /// m-squared
      TH2D *mass2De = new TH2D("mass2De", "", nbins, 0., 4., nbins, -0.5, 10.);
      fMSquaredMap.insert(pair<MpdPidUtils::ePartType, TH2D *>(MpdPidUtils::kDeuteron, mass2De));

      /// Abundance
      vecTH1Dptrs vecYield;
      TH1D       *hYield_de = new TH1D("hYield_de", "", 300, 0., 3.);
      vecYield.push_back(hYield_de);
      fAbundanceMap.insert(pair<MpdPidUtils::ePartType, vecTH1Dptrs>(MpdPidUtils::kDeuteron, vecYield));

      /// Efficiency & Contamination
      TH1D                    *hEC[4];
      vecTH1Dptrs              vecEC;
      std::vector<vecTH1Dptrs> vecvecEC;
      for (Int_t iHist = 0; iHist < 4; iHist++) {
         hEC[iHist] = new TH1D(ECHistNames[iHist] + "_" + MpdPidUtils::cParticleShortName[MpdPidUtils::kDeuteron], "",
                               50, PMIN, PMAX);
         hEC[iHist]->Sumw2();
         vecEC.push_back(hEC[iHist]);
      }
      vecvecEC.push_back(vecEC);
      fEffContMap.insert(pair<MpdPidUtils::ePartType, std::vector<vecTH1Dptrs>>(MpdPidUtils::kDeuteron, vecvecEC));
   }
   if (Particles.Contains("tr")) {
      /// dE/dx
      if (fTrackingState == MpdPidUtils::kCFHM) {
         MaxEnLossLoc = 120.;
         for (Int_t i = 0; i < 4; i++) XMIN[i][MpdPidUtils::kTriton] = 20.;
         for (Int_t i = 4; i < 8; i++) XMIN[i][MpdPidUtils::kTriton] = 15.;
         XMIN[8][MpdPidUtils::kTriton]  = 10.;
         XMIN[9][MpdPidUtils::kTriton]  = 10.;
         XMIN[10][MpdPidUtils::kTriton] = 8.;
         XMIN[11][MpdPidUtils::kTriton] = 7.;
         XMIN[12][MpdPidUtils::kTriton] = 6.;
         XMIN[13][MpdPidUtils::kTriton] = 6.;
         XMIN[14][MpdPidUtils::kTriton] = 5.;
         XMIN[15][MpdPidUtils::kTriton] = 5.;
         XMIN[16][MpdPidUtils::kTriton] = 5.;
         XMIN[17][MpdPidUtils::kTriton] = 5.;
         XMIN[18][MpdPidUtils::kTriton] = 4.;
         XMIN[19][MpdPidUtils::kTriton] = 4.;
         XMIN[20][MpdPidUtils::kTriton] = 4.;
         XMIN[21][MpdPidUtils::kTriton] = 3.;
         XMIN[22][MpdPidUtils::kTriton] = 3.;
         XMIN[23][MpdPidUtils::kTriton] = 3.;
         XMIN[24][MpdPidUtils::kTriton] = 2.5;
         XMIN[25][MpdPidUtils::kTriton] = 2.5;
         XMIN[26][MpdPidUtils::kTriton] = 2.5;
         XMIN[27][MpdPidUtils::kTriton] = 2.5;
         for (Int_t i = 28; i < 33; i++) XMIN[i][MpdPidUtils::kTriton] = 2.;
         for (Int_t i = 33; i < 37; i++) XMIN[i][MpdPidUtils::kTriton] = 1.5;
         XMAX[0][MpdPidUtils::kTriton]  = 120.;
         XMAX[1][MpdPidUtils::kTriton]  = 120.;
         XMAX[2][MpdPidUtils::kTriton]  = 100.;
         XMAX[3][MpdPidUtils::kTriton]  = 80.;
         XMAX[4][MpdPidUtils::kTriton]  = 70.;
         XMAX[5][MpdPidUtils::kTriton]  = 60.;
         XMAX[6][MpdPidUtils::kTriton]  = 50.;
         XMAX[7][MpdPidUtils::kTriton]  = 45.;
         XMAX[8][MpdPidUtils::kTriton]  = 40.;
         XMAX[9][MpdPidUtils::kTriton]  = 35.;
         XMAX[10][MpdPidUtils::kTriton] = 30.;
         XMAX[11][MpdPidUtils::kTriton] = 26.;
         XMAX[12][MpdPidUtils::kTriton] = 24.;
         XMAX[13][MpdPidUtils::kTriton] = 21.;
         XMAX[14][MpdPidUtils::kTriton] = 19.;
         XMAX[15][MpdPidUtils::kTriton] = 17.;
         XMAX[16][MpdPidUtils::kTriton] = 15.;
         XMAX[17][MpdPidUtils::kTriton] = 13.;
         XMAX[18][MpdPidUtils::kTriton] = 13.;
         XMAX[19][MpdPidUtils::kTriton] = 12.;
         XMAX[20][MpdPidUtils::kTriton] = 11.;
         XMAX[21][MpdPidUtils::kTriton] = 11.;
         XMAX[22][MpdPidUtils::kTriton] = 10.;
         XMAX[23][MpdPidUtils::kTriton] = 9.;
         XMAX[24][MpdPidUtils::kTriton] = 8.5;
         XMAX[25][MpdPidUtils::kTriton] = 8.;
         XMAX[26][MpdPidUtils::kTriton] = 7.5;
         XMAX[27][MpdPidUtils::kTriton] = 7.;
         XMAX[28][MpdPidUtils::kTriton] = 7.;
         XMAX[29][MpdPidUtils::kTriton] = 6.5;
         XMAX[30][MpdPidUtils::kTriton] = 6.;
         XMAX[31][MpdPidUtils::kTriton] = 5.5;
         XMAX[32][MpdPidUtils::kTriton] = 5.;
         XMAX[33][MpdPidUtils::kTriton] = 4.5;
         XMAX[34][MpdPidUtils::kTriton] = 4.5;
         XMAX[35][MpdPidUtils::kTriton] = 4.;
         XMAX[36][MpdPidUtils::kTriton] = 4.;
      } else if (fTrackingState == MpdPidUtils::kCF) {
         MaxEnLossLoc                  = 500000.;
         XMIN[0][MpdPidUtils::kTriton] = 100000.;
         XMIN[1][MpdPidUtils::kTriton] = 50000.;
         XMIN[2][MpdPidUtils::kTriton] = 40000.;
         XMIN[3][MpdPidUtils::kTriton] = 30000.;
         XMIN[4][MpdPidUtils::kTriton] = 20000.;
         XMIN[5][MpdPidUtils::kTriton] = 10000.;
         XMIN[6][MpdPidUtils::kTriton] = 8000.;
         XMIN[7][MpdPidUtils::kTriton] = 5000.;
         for (Int_t i = 8; i < 14; i++) XMIN[i][MpdPidUtils::kTriton] = 2000.;
         for (Int_t i = 14; i < 20; i++) XMIN[i][MpdPidUtils::kTriton] = 5000.;
         for (Int_t i = 20; i < 26; i++) XMIN[i][MpdPidUtils::kTriton] = 6000.;
         for (Int_t i = 26; i < 29; i++) XMIN[i][MpdPidUtils::kTriton] = 5000.;
         for (Int_t i = 29; i < 35; i++) XMIN[i][MpdPidUtils::kTriton] = 4000.;
         XMIN[35][MpdPidUtils::kTriton] = 3500.;
         XMIN[36][MpdPidUtils::kTriton] = 3500.;
         XMAX[0][MpdPidUtils::kTriton]  = 500000.;
         XMAX[1][MpdPidUtils::kTriton]  = 500000.;
         XMAX[2][MpdPidUtils::kTriton]  = 350000.;
         XMAX[3][MpdPidUtils::kTriton]  = 250000.;
         XMAX[4][MpdPidUtils::kTriton]  = 200000.;
         XMAX[5][MpdPidUtils::kTriton]  = 160000.;
         XMAX[6][MpdPidUtils::kTriton]  = 140000.;
         XMAX[7][MpdPidUtils::kTriton]  = 100000.;
         XMAX[8][MpdPidUtils::kTriton]  = 90000.;
         XMAX[9][MpdPidUtils::kTriton]  = 80000.;
         XMAX[10][MpdPidUtils::kTriton] = 60000.;
         XMAX[11][MpdPidUtils::kTriton] = 50000.;
         XMAX[12][MpdPidUtils::kTriton] = 45000.;
         XMAX[13][MpdPidUtils::kTriton] = 40000.;
         XMAX[14][MpdPidUtils::kTriton] = 40000.;
         XMAX[15][MpdPidUtils::kTriton] = 30000.;
         XMAX[16][MpdPidUtils::kTriton] = 30000.;
         XMAX[17][MpdPidUtils::kTriton] = 30000.;
         XMAX[18][MpdPidUtils::kTriton] = 25000.;
         XMAX[19][MpdPidUtils::kTriton] = 25000.;
         XMAX[20][MpdPidUtils::kTriton] = 22000.;
         XMAX[21][MpdPidUtils::kTriton] = 22000.;
         XMAX[22][MpdPidUtils::kTriton] = 20000.;
         XMAX[23][MpdPidUtils::kTriton] = 18000.;
         XMAX[24][MpdPidUtils::kTriton] = 18000.;
         XMAX[25][MpdPidUtils::kTriton] = 16000.;
         XMAX[26][MpdPidUtils::kTriton] = 16000.;
         XMAX[27][MpdPidUtils::kTriton] = 15000.;
         XMAX[28][MpdPidUtils::kTriton] = 14000.;
         XMAX[29][MpdPidUtils::kTriton] = 14000.;
         XMAX[30][MpdPidUtils::kTriton] = 13000.;
         XMAX[31][MpdPidUtils::kTriton] = 11000.;
         XMAX[32][MpdPidUtils::kTriton] = 11000.;
         XMAX[33][MpdPidUtils::kTriton] = 10000.;
         XMAX[34][MpdPidUtils::kTriton] = 10000.;
         XMAX[35][MpdPidUtils::kTriton] = 9000.;
         XMAX[36][MpdPidUtils::kTriton] = 8000.;
      } else {
         MaxEnLossLoc                  = 200.e-06;
         XMIN[0][MpdPidUtils::kTriton] = 3.e-06;
         XMIN[1][MpdPidUtils::kTriton] = 3.e-06;
         XMIN[2][MpdPidUtils::kTriton] = 2.e-06;
         for (Int_t i = 0; i < (MpdPidUtils::nQAHists - 3); i++) XMIN[i][MpdPidUtils::kTriton] = 1.e-06;
         XMAX[0][MpdPidUtils::kTriton] = 200.e-06;
         XMAX[1][MpdPidUtils::kTriton] = 160.e-06;
         XMAX[2][MpdPidUtils::kTriton] = 140.e-06;
         XMAX[3][MpdPidUtils::kTriton] = 110.e-06;
         XMAX[4][MpdPidUtils::kTriton] = 80.e-06;
         XMAX[5][MpdPidUtils::kTriton] = 60.e-06;
         XMAX[6][MpdPidUtils::kTriton] = 40.e-06;
         for (Int_t i = 7; i < (MpdPidUtils::nQAHists - 3); i++) XMAX[i][MpdPidUtils::kTriton] = 20.e-06;
      }
      MpdPidUtils::dEdXStruct_t STr;
      for (Int_t iHist = 0; iHist < MpdPidUtils::nQAHists; iHist++) {
         histName = "Tr_" + TString::Itoa(iHist, 10);
         STr.dEdXPart[iHist] =
            new TH1D(histName, "", nbins, XMIN[iHist][MpdPidUtils::kTriton], XMAX[iHist][MpdPidUtils::kTriton]);
      }
      STr.BetheBlochHist = new TH2D("hTrBB_QA", "", 300, 0.05, 3.0, 300, MinEnLoss, MaxEnLossLoc);
      STr.ibeg           = 3;
      STr.iend           = 39;
      fEnLossMap.insert(pair<MpdPidUtils::ePartType, MpdPidUtils::dEdXStruct_t>(MpdPidUtils::kTriton, STr));

      /// m-squared
      TH2D *mass2Tr = new TH2D("mass2Tr", "", nbins, 0., 4., nbins, -0.5, 15.);
      fMSquaredMap.insert(pair<MpdPidUtils::ePartType, TH2D *>(MpdPidUtils::kTriton, mass2Tr));

      /// Abundance
      vecTH1Dptrs vecYield;
      TH1D       *hYield_tr = new TH1D("hYield_tr", "", 300, 0., 3.);
      vecYield.push_back(hYield_tr);
      fAbundanceMap.insert(pair<MpdPidUtils::ePartType, vecTH1Dptrs>(MpdPidUtils::kTriton, vecYield));

      /// Efficiency & Contamination
      TH1D                    *hEC[4];
      vecTH1Dptrs              vecEC;
      std::vector<vecTH1Dptrs> vecvecEC;
      for (Int_t iHist = 0; iHist < 4; iHist++) {
         hEC[iHist] = new TH1D(ECHistNames[iHist] + "_" + MpdPidUtils::cParticleShortName[MpdPidUtils::kTriton], "", 50,
                               PMIN, PMAX);
         hEC[iHist]->Sumw2();
         vecEC.push_back(hEC[iHist]);
      }
      vecvecEC.push_back(vecEC);
      fEffContMap.insert(pair<MpdPidUtils::ePartType, std::vector<vecTH1Dptrs>>(MpdPidUtils::kTriton, vecvecEC));
   }
   if (Particles.Contains("he3")) {
      /// dE/dx
      if (fTrackingState == MpdPidUtils::kCFHM) {
         MaxEnLossLoc = 150.;
         for (Int_t i = 0; i < 5; i++) XMIN[i][MpdPidUtils::kHe3] = 5.;
         XMIN[5][MpdPidUtils::kHe3] = 10.;
         XMIN[6][MpdPidUtils::kHe3] = 10.;
         XMIN[7][MpdPidUtils::kHe3] = 15.;
         XMIN[8][MpdPidUtils::kHe3] = 15.;
         for (Int_t i = 9; i < 15; i++) XMIN[i][MpdPidUtils::kHe3] = 10.;
         for (Int_t i = 15; i < 19; i++) XMIN[i][MpdPidUtils::kHe3] = 8.;
         for (Int_t i = 19; i < 26; i++) XMIN[i][MpdPidUtils::kHe3] = 7.;
         XMIN[26][MpdPidUtils::kHe3] = 6.;
         XMIN[27][MpdPidUtils::kHe3] = 6.;
         XMAX[0][MpdPidUtils::kHe3]  = 150.;
         XMAX[1][MpdPidUtils::kHe3]  = 150.;
         XMAX[2][MpdPidUtils::kHe3]  = 150.;
         XMAX[3][MpdPidUtils::kHe3]  = 120.;
         XMAX[4][MpdPidUtils::kHe3]  = 100.;
         XMAX[5][MpdPidUtils::kHe3]  = 80.;
         XMAX[6][MpdPidUtils::kHe3]  = 70.;
         XMAX[7][MpdPidUtils::kHe3]  = 60.;
         XMAX[8][MpdPidUtils::kHe3]  = 55.;
         XMAX[9][MpdPidUtils::kHe3]  = 50.;
         XMAX[10][MpdPidUtils::kHe3] = 45.;
         XMAX[11][MpdPidUtils::kHe3] = 40.;
         XMAX[12][MpdPidUtils::kHe3] = 35.;
         XMAX[13][MpdPidUtils::kHe3] = 30.;
         XMAX[14][MpdPidUtils::kHe3] = 30.;
         XMAX[15][MpdPidUtils::kHe3] = 25.;
         XMAX[16][MpdPidUtils::kHe3] = 23.;
         XMAX[17][MpdPidUtils::kHe3] = 21.;
         XMAX[18][MpdPidUtils::kHe3] = 20.;
         XMAX[19][MpdPidUtils::kHe3] = 19.;
         XMAX[20][MpdPidUtils::kHe3] = 18.;
         XMAX[21][MpdPidUtils::kHe3] = 17.;
         XMAX[22][MpdPidUtils::kHe3] = 16.;
         XMAX[23][MpdPidUtils::kHe3] = 16.;
         XMAX[24][MpdPidUtils::kHe3] = 16.;
         XMAX[25][MpdPidUtils::kHe3] = 15.;
         XMAX[26][MpdPidUtils::kHe3] = 15.;
         XMAX[27][MpdPidUtils::kHe3] = 15.;
      } else if (fTrackingState == MpdPidUtils::kCF) {
         MaxEnLossLoc               = 1800000.;
         XMIN[0][MpdPidUtils::kHe3] = 100000.;
         XMIN[1][MpdPidUtils::kHe3] = 100000.;
         XMIN[2][MpdPidUtils::kHe3] = 50000.;
         XMIN[3][MpdPidUtils::kHe3] = 30000.;
         for (Int_t i = 4; i < (MpdPidUtils::nQAHists - 9); i++) XMIN[i][MpdPidUtils::kHe3] = 10000.;
         XMAX[0][MpdPidUtils::kHe3]  = 1800000.;
         XMAX[1][MpdPidUtils::kHe3]  = 1000000.;
         XMAX[2][MpdPidUtils::kHe3]  = 800000.;
         XMAX[3][MpdPidUtils::kHe3]  = 500000.;
         XMAX[4][MpdPidUtils::kHe3]  = 350000.;
         XMAX[5][MpdPidUtils::kHe3]  = 250000.;
         XMAX[6][MpdPidUtils::kHe3]  = 200000.;
         XMAX[7][MpdPidUtils::kHe3]  = 180000.;
         XMAX[8][MpdPidUtils::kHe3]  = 150000.;
         XMAX[9][MpdPidUtils::kHe3]  = 130000.;
         XMAX[10][MpdPidUtils::kHe3] = 110000.;
         XMAX[11][MpdPidUtils::kHe3] = 90000.;
         XMAX[12][MpdPidUtils::kHe3] = 80000.;
         XMAX[13][MpdPidUtils::kHe3] = 70000.;
         XMAX[14][MpdPidUtils::kHe3] = 60000.;
         XMAX[15][MpdPidUtils::kHe3] = 50000.;
         XMAX[16][MpdPidUtils::kHe3] = 45000.;
         XMAX[17][MpdPidUtils::kHe3] = 45000.;
         XMAX[18][MpdPidUtils::kHe3] = 40000.;
         XMAX[19][MpdPidUtils::kHe3] = 40000.;
         XMAX[20][MpdPidUtils::kHe3] = 35000.;
         XMAX[21][MpdPidUtils::kHe3] = 35000.;
         XMAX[22][MpdPidUtils::kHe3] = 30000.;
         XMAX[23][MpdPidUtils::kHe3] = 30000.;
         XMAX[24][MpdPidUtils::kHe3] = 28000.;
         XMAX[25][MpdPidUtils::kHe3] = 28000.;
         XMAX[26][MpdPidUtils::kHe3] = 26000.;
         XMAX[27][MpdPidUtils::kHe3] = 26000.;
         XMAX[28][MpdPidUtils::kHe3] = 30000.;
         XMAX[29][MpdPidUtils::kHe3] = 30000.;
         XMAX[30][MpdPidUtils::kHe3] = 30000.;
      } else {
         MaxEnLossLoc               = 450.e-06;
         XMIN[0][MpdPidUtils::kHe3] = 3.e-06;
         XMIN[1][MpdPidUtils::kHe3] = 3.e-06;
         XMIN[2][MpdPidUtils::kHe3] = 3.e-06;
         XMIN[3][MpdPidUtils::kHe3] = 2.e-06;
         XMIN[4][MpdPidUtils::kHe3] = 2.e-06;
         for (Int_t i = 5; i < 9; i++) XMIN[i][MpdPidUtils::kHe3] = 1.5e-06;
         for (Int_t i = 9; i < (MpdPidUtils::nQAHists - 9); i++) XMIN[i][MpdPidUtils::kHe3] = 1.e-06;
         XMAX[0][MpdPidUtils::kHe3]  = 450.e-06;
         XMAX[1][MpdPidUtils::kHe3]  = 350.e-06;
         XMAX[2][MpdPidUtils::kHe3]  = 300.e-06;
         XMAX[3][MpdPidUtils::kHe3]  = 250.e-06;
         XMAX[4][MpdPidUtils::kHe3]  = 200.e-06;
         XMAX[5][MpdPidUtils::kHe3]  = 150.e-06;
         XMAX[6][MpdPidUtils::kHe3]  = 100.e-06;
         XMAX[7][MpdPidUtils::kHe3]  = 80.e-06;
         XMAX[8][MpdPidUtils::kHe3]  = 70.e-06;
         XMAX[9][MpdPidUtils::kHe3]  = 60.e-06;
         XMAX[10][MpdPidUtils::kHe3] = 50.e-06;
         XMAX[11][MpdPidUtils::kHe3] = 40.e-06;
         XMAX[12][MpdPidUtils::kHe3] = 35.e-06;
         for (Int_t i = 13; i < (MpdPidUtils::nQAHists - 9); i++)
            XMAX[i][MpdPidUtils::kHe3] = XMAX[11][MpdPidUtils::kHe3] = 32.e-06;
      }
      MpdPidUtils::dEdXStruct_t SHe3;
      for (Int_t iHist = 0; iHist < MpdPidUtils::nQAHists; iHist++) {
         histName = "He3_" + TString::Itoa(iHist, 10);
         SHe3.dEdXPart[iHist] =
            new TH1D(histName, "", nbins, XMIN[iHist][MpdPidUtils::kHe3], XMAX[iHist][MpdPidUtils::kHe3]);
      }
      SHe3.BetheBlochHist = new TH2D("hHe3BB_QA", "", 300, 0.05, 3.0, 300, MinEnLoss, MaxEnLossLoc);
      SHe3.ibeg           = (fTrackingState == MpdPidUtils::kCFHM) ? 3 : 2;
      SHe3.iend           = (fTrackingState == MpdPidUtils::kCFHM) ? 30 : 32;
      fEnLossMap.insert(pair<MpdPidUtils::ePartType, MpdPidUtils::dEdXStruct_t>(MpdPidUtils::kHe3, SHe3));

      /// m-squared
      TH2D *mass2He3 = new TH2D("mass2He3", "", nbins, 0., 4., nbins, -0.5, 12.);
      fMSquaredMap.insert(pair<MpdPidUtils::ePartType, TH2D *>(MpdPidUtils::kHe3, mass2He3));

      /// Abundance
      vecTH1Dptrs vecYield;
      TH1D       *hYield_he3 = new TH1D("hYield_he3", "", 150, 0., 1.5);
      vecYield.push_back(hYield_he3);
      fAbundanceMap.insert(pair<MpdPidUtils::ePartType, vecTH1Dptrs>(MpdPidUtils::kHe3, vecYield));

      /// Efficiency & Contamination
      TH1D                    *hEC[4];
      vecTH1Dptrs              vecEC;
      std::vector<vecTH1Dptrs> vecvecEC;
      for (Int_t iHist = 0; iHist < 4; iHist++) {
         hEC[iHist] =
            new TH1D(ECHistNames[iHist] + "_" + MpdPidUtils::cParticleShortName[MpdPidUtils::kHe3], "", 50, PMIN, PMAX);
         hEC[iHist]->Sumw2();
         vecEC.push_back(hEC[iHist]);
      }
      vecvecEC.push_back(vecEC);
      fEffContMap.insert(pair<MpdPidUtils::ePartType, std::vector<vecTH1Dptrs>>(MpdPidUtils::kHe3, vecvecEC));
   }
   if (Particles.Contains("he4")) {
      /// dE/dx
      if (fTrackingState == MpdPidUtils::kCFHM) {
         MaxEnLossLoc = 250.;
         for (Int_t i = 0; i < 4; i++) XMIN[i][8] = 5.;
         XMIN[4][8]  = 20.;
         XMIN[5][8]  = 20.;
         XMIN[6][8]  = 25.;
         XMIN[7][8]  = 25.;
         XMIN[8][8]  = 25.;
         XMIN[9][8]  = 20.;
         XMIN[10][8] = 20.;
         XMIN[11][8] = 20.;
         XMIN[12][8] = 18.;
         XMIN[13][8] = 15.;
         XMIN[14][8] = 15.;
         XMIN[15][8] = 14.;
         XMIN[16][8] = 14.;
         XMIN[17][8] = 13.;
         XMIN[18][8] = 12.;
         XMIN[19][8] = 11.;
         XMIN[26][8] = 9.;
         XMIN[27][8] = 8.;
         for (Int_t i = 20; i < 26; i++) XMIN[i][8] = 10.;
         XMAX[0][8]  = 250.;
         XMAX[1][8]  = 200.;
         XMAX[2][8]  = 150.;
         XMAX[3][8]  = 140.;
         XMAX[4][8]  = 120.;
         XMAX[5][8]  = 100.;
         XMAX[6][8]  = 90.;
         XMAX[7][8]  = 80.;
         XMAX[8][8]  = 75.;
         XMAX[9][8]  = 70.;
         XMAX[10][8] = 60.;
         XMAX[11][8] = 55.;
         XMAX[12][8] = 50.;
         XMAX[13][8] = 45.;
         XMAX[14][8] = 42.;
         XMAX[15][8] = 40.;
         XMAX[16][8] = 35.;
         XMAX[17][8] = 33.;
         XMAX[18][8] = 30.;
         XMAX[19][8] = 28.;
         XMAX[20][8] = 26.;
         XMAX[21][8] = 25.;
         XMAX[22][8] = 23.;
         XMAX[23][8] = 21.;
         XMAX[24][8] = 20.;
         XMAX[25][8] = 19.;
         XMAX[26][8] = 19.;
         XMAX[27][8] = 22.;
      } else if (fTrackingState == MpdPidUtils::kCF) {
         MaxEnLossLoc = 1800000.;
         XMIN[0][8]   = 100000.;
         XMIN[1][8]   = 80000.;
         XMIN[2][8]   = 50000.;
         XMIN[3][8]   = 30000.;
         XMIN[4][8]   = 30000.;
         XMIN[5][8]   = 20000.;
         for (Int_t i = 6; i < (MpdPidUtils::nQAHists - 12); i++) XMIN[i][8] = 10000.;
         XMAX[0][8]  = 1800000.;
         XMAX[1][8]  = 1200000.;
         XMAX[2][8]  = 900000.;
         XMAX[3][8]  = 700000.;
         XMAX[4][8]  = 500000.;
         XMAX[5][8]  = 350000.;
         XMAX[6][8]  = 300000.;
         XMAX[7][8]  = 250000.;
         XMAX[8][8]  = 200000.;
         XMAX[9][8]  = 180000.;
         XMAX[10][8] = 160000.;
         XMAX[11][8] = 140000.;
         XMAX[12][8] = 130000.;
         XMAX[13][8] = 120000.;
         XMAX[14][8] = 100000.;
         XMAX[15][8] = 90000.;
         XMAX[16][8] = 80000.;
         XMAX[17][8] = 70000.;
         XMAX[18][8] = 70000.;
         XMAX[19][8] = 60000.;
         XMAX[20][8] = 55000.;
         XMAX[21][8] = 50000.;
         XMAX[22][8] = 45000.;
         XMAX[23][8] = 45000.;
         XMAX[24][8] = 40000.;
         XMAX[25][8] = 40000.;
         XMAX[26][8] = 40000.;
         XMAX[27][8] = 40000.;
      } else {
         MaxEnLossLoc = 350.e-06;
         XMIN[0][8]   = 3.e-06;
         XMIN[1][8]   = 2.e-06;
         for (Int_t i = 2; i < (MpdPidUtils::nQAHists - 12); i++) XMIN[i][8] = 1.e-06;
         XMAX[0][8]  = 350.e-06;
         XMAX[1][8]  = 300.e-06;
         XMAX[2][8]  = 250.e-06;
         XMAX[3][8]  = 200.e-06;
         XMAX[4][8]  = 150.e-06;
         XMAX[5][8]  = 100.e-06;
         XMAX[6][8]  = 80.e-06;
         XMAX[7][8]  = 70.e-06;
         XMAX[8][8]  = 60.e-06;
         XMAX[9][8]  = 50.e-06;
         XMAX[10][8] = 40.e-06;
         for (Int_t i = 11; i < (MpdPidUtils::nQAHists - 12); i++) XMAX[i][8] = 32.e-06;
      }
      MpdPidUtils::dEdXStruct_t SHe4;
      for (Int_t iHist = 0; iHist < MpdPidUtils::nQAHists; iHist++) {
         histName             = "He4_" + TString::Itoa(iHist, 10);
         SHe4.dEdXPart[iHist] = new TH1D(histName, "", (iHist < 27) ? nbins : 50, XMIN[iHist][8], XMAX[iHist][8]);
      }
      SHe4.BetheBlochHist = new TH2D("hHe4BB_QA", "", 300, 0.05, 3.0, 300, MinEnLoss, MaxEnLossLoc);
      SHe4.ibeg           = 3;
      SHe4.iend           = 30;
      fEnLossMap.insert(pair<MpdPidUtils::ePartType, MpdPidUtils::dEdXStruct_t>(MpdPidUtils::kHe4, SHe4));

      /// m-squared
      TH2D *mass2He4 = new TH2D("mass2He4", "", nbins, 0., 4., nbins, -0.5, 12.);
      fMSquaredMap.insert(pair<MpdPidUtils::ePartType, TH2D *>(MpdPidUtils::kHe4, mass2He4));

      /// Abundance
      vecTH1Dptrs vecYield;
      TH1D       *hYield_he4 = new TH1D("hYield_he4", "", 150, 0., 1.5);
      vecYield.push_back(hYield_he4);
      fAbundanceMap.insert(pair<MpdPidUtils::ePartType, vecTH1Dptrs>(MpdPidUtils::kHe4, vecYield));

      /// Efficiency & Contamination
      TH1D                    *hEC[4];
      vecTH1Dptrs              vecEC;
      std::vector<vecTH1Dptrs> vecvecEC;
      for (Int_t iHist = 0; iHist < 4; iHist++) {
         hEC[iHist] =
            new TH1D(ECHistNames[iHist] + "_" + MpdPidUtils::cParticleShortName[MpdPidUtils::kHe4], "", 50, PMIN, PMAX);
         hEC[iHist]->Sumw2();
         vecEC.push_back(hEC[iHist]);
      }
      vecvecEC.push_back(vecEC);
      fEffContMap.insert(pair<MpdPidUtils::ePartType, std::vector<vecTH1Dptrs>>(MpdPidUtils::kHe4, vecvecEC));
   }

   if (fTrackingState == MpdPidUtils::kCFHM)
      MaxEnLoss = 60.;
   else if (fTrackingState == MpdPidUtils::kCF)
      MaxEnLoss = 60000.;
   else
      MaxEnLoss = 60.e-06;

   fSumBetheBlochHists = new TH2D("fSumBetheBlochHists", "", 300, 0.05, 3.0, 300, MinEnLoss, MaxEnLoss);
   fChBetheBlochHists  = new TH2D("fChBetheBlochHists", "", 600, -1.5, 1.5, 600, MinEnLoss, MaxEnLoss);
   fm2LightHist        = new TH2D("fm2LightHist", "", 200, 0.0, 3.0, 200, -0.5, 1.5);
   fm2HeavyHist        = new TH2D("fm2HeavyHist", "", 200, 0.0, 5.0, 200, -0.5, 12.0);
   for (Int_t i = 0; i < (MpdPidUtils::nQAHists - 8); i++) Xlow[i] = 0.05 * (1. + i);
   Xlow[32] = 1.7;
   Xlow[33] = 1.8;
   Xlow[34] = 1.9;
   Xlow[35] = 2.0;
   Xlow[36] = 2.1;
   Xlow[37] = 2.2;
   Xlow[38] = 2.4;
   Xlow[39] = 2.6;
   for (Int_t i = 0; i < (MpdPidUtils::nQAHists - 9); i++) Xhigh[i] = 0.05 + 0.05 * (1. + i);
   Xhigh[31] = 1.7;
   Xhigh[32] = 1.8;
   Xhigh[33] = 1.9;
   Xhigh[34] = 2.0;
   Xhigh[35] = 2.1;
   Xhigh[36] = 2.2;
   Xhigh[37] = 2.4;
   Xhigh[38] = 2.6;
   Xhigh[39] = 3.0;
   for (Int_t i = 0; i < 40; i++) X[i] = (Xlow[i] + Xhigh[i]) / 2;
}

MpdPidUtils::ePartType MpdPidQA::GetPartType(Int_t pdg)
{
   map<Int_t, MpdPidUtils::ePartType>::iterator ret = fPartTypeMap.find(pdg);
   if (ret != fPartTypeMap.end())
      return ret->second;
   else
      return MpdPidUtils::kUnknown;
}

void MpdPidQA::FillDedxHists(Double_t p, Double_t dedx, Int_t pdg)
{
   Double_t sign = pdg > 0 ? 1.0 : -1.0;
   auto     ret  = fEnLossMap.find(GetPartType(pdg));
   if (ret != fEnLossMap.end()) {
      fSumBetheBlochHists->Fill(p, dedx);
      fChBetheBlochHists->Fill(sign * p, dedx);
      ret->second.BetheBlochHist->Fill(p, dedx);
      for (Int_t i = 0; i < MpdPidUtils::nQAHists; i++) {
         if ((p > Xlow[i]) && (p <= Xhigh[i])) {
            if ((i < ret->second.ibeg) || (i > ret->second.iend)) continue;
            ret->second.dEdXPart[i - ret->second.ibeg]->Fill(dedx);
         }
      }
   }
}

void MpdPidQA::Fillm2Hists(Double_t p, Double_t m2, Int_t pdg)
{
   Int_t pdgc = TMath::Abs(pdg);
   auto  ret  = fMSquaredMap.find(GetPartType(pdg));
   if (ret != fMSquaredMap.end()) ret->second->Fill(p, m2);

   if ((pdgc == 211) || (pdgc == 321) || (pdgc == 2212)) fm2LightHist->Fill(p, m2);
   if ((pdgc == PDG_DEUTERON) || (pdgc == PDG_TRITON) || (pdgc == PDG_HE3) || (pdgc == PDG_HE4))
      fm2HeavyHist->Fill(p, m2);
}

void MpdPidQA::FillAmplHists(Double_t p, Int_t pdg)
{
   if (p <= 5.0) {
      MpdPidUtils::ePartCharge ech;
      if ((TMath::Abs(pdg) == 11) || (TMath::Abs(pdg) == 13))
         ech = pdg > 0 ? MpdPidUtils::kNeg : MpdPidUtils::kPos;
      else
         ech = pdg > 0 ? MpdPidUtils::kPos : MpdPidUtils::kNeg;
      auto ret = fAbundanceMap.find(GetPartType(pdg));
      if (ret != fAbundanceMap.end()) ret->second[ech]->Fill(p);
   }
}

Bool_t MpdPidQA::FillEffContHists(Double_t p, Double_t dedx, Int_t charge, Int_t pdg, Double_t fProbCut)
{
   FillEffDenominator(p, pdg);
   if (!FillProbs(p, dedx, charge)) return kFALSE;

   Int_t    pidpdg  = GetMaxProb();
   Double_t maxprob = GetProb(GetPartType(pidpdg));

   if (maxprob < fProbCut) return kFALSE;
   MpdPidUtils::ePartCharge ech = charge > 0 ? MpdPidUtils::kPos : MpdPidUtils::kNeg;
   FillEffContHists(p, pidpdg, pdg, ech);
   return kTRUE;
}

Bool_t MpdPidQA::FillEffContHists(Double_t p, Double_t dedx, Double_t m2, Int_t charge, Int_t pdg, Double_t fProbCut)
{
   FillEffDenominator(p, pdg);
   if (!FillProbs(p, dedx, m2, charge)) return kFALSE;

   Int_t    pidpdg  = GetMaxProb();
   Double_t maxprob = GetProb(GetPartType(pidpdg));

   if (maxprob < fProbCut) return kFALSE;
   MpdPidUtils::ePartCharge ech = charge > 0 ? MpdPidUtils::kPos : MpdPidUtils::kNeg;
   FillEffContHists(p, pidpdg, pdg, ech);
   return kTRUE;
}

Bool_t MpdPidQA::FillEffContHists(Double_t p, Int_t pidpdg, Int_t pdg, MpdPidUtils::ePartCharge ech)
{
   auto ret = fEffContMap.find(GetPartType(pidpdg));
   if (ret != fEffContMap.end()) {
      ret->second[ech][1]->Fill(p);
      if (pidpdg == pdg)
         ret->second[ech][2]->Fill(p);
      else
         ret->second[ech][3]->Fill(p);
   }
   return kTRUE;
}

void MpdPidQA::FillEffDenominator(Double_t p, Int_t pdg)
{
   auto                     ret = fEffContMap.find(GetPartType(pdg));
   MpdPidUtils::ePartCharge ech;
   if ((TMath::Abs(pdg) == 11) || (TMath::Abs(pdg) == 13))
      ech = pdg > 0 ? MpdPidUtils::kNeg : MpdPidUtils::kPos;
   else
      ech = pdg > 0 ? MpdPidUtils::kPos : MpdPidUtils::kNeg;
   if (ret != fEffContMap.end()) ret->second[ech][0]->Fill(p);
}

void MpdPidQA::GetDedxQA(TString dir)
{
   Double_t                                          XFUNCMIN, XFUNCMAX, mom, intF_MultX, intF;
   TString                                           FName, FMXName;
   Int_t                                             ibegloc, iendloc, nQAHistsloc, inam = 0;
   map<MpdPidUtils::ePartType, TGraphAsymmErrors *>  graphs;
   map<MpdPidUtils::ePartType, vecTF1ptrs>::iterator it_parBB;
   vecTF1ptrs                                        vecParBB, vecParBBMultX;
   /*
   for ( auto it = fEnLossMap.begin(); it != fEnLossMap.end(); ++it ) {
      ibegloc = it->second.ibeg; iendloc = it->second.iend; nQAHistsloc = iendloc - ibegloc + 1;
      it_parBB = fdEdxBBMap.find( it->first );
      vecParBB = it_parBB->second;

      for ( TF1* BBF : vecParBB ) {
         FName = BBF->GetName(); FName.Append("_MultX");
         TF1* parBBMultX = new TF1( FName,[&](Double_t *x, Double_t *p){ return p[0] * x[0] * BBF->Eval(x[0]);
   },BBF->GetXmin(),BBF->GetXmax(),1); parBBMultX->SetParameter(0, 1.0); vecParBBMultX.push_back(parBBMultX);
      }

      Double_t *Xgraph = new Double_t[nQAHistsloc]; Double_t *Xerrgraph = new Double_t[nQAHistsloc];
      Double_t *Ygraph = new Double_t[nQAHistsloc]; Double_t *sigma = new Double_t[nQAHistsloc];
      Double_t *sigma1 = new Double_t[nQAHistsloc]; Double_t *sigma2 = new Double_t[nQAHistsloc];
      for ( Int_t k = 0; k < nQAHistsloc; k++ ) { Xgraph[k] = X[k + ibegloc]; Xerrgraph[k] = 0.; }

      for ( Int_t i = 0; i < nQAHistsloc; i++ ) {
         XFUNCMIN = it->second.dEdXPart[i]->GetXaxis()->GetXmin(), XFUNCMAX =
   it->second.dEdXPart[i]->GetXaxis()->GetXmax(); TF1* Gaus = new TF1("Gaus", "gaus", XFUNCMIN, XFUNCMAX);
         it->second.dEdXPart[i]->Fit(Gaus, "Q0R");
         TF1* AGaus = new TF1("AGaus", this, &MpdPid::AsymGaus, XFUNCMIN, XFUNCMAX, 4, "MpdPid", "AsymGaus");
         AGaus->SetParameters(Gaus->GetParameter(0), Gaus->GetParameter(1), Gaus->GetParameter(2), 0.01);
         it->second.dEdXPart[i]->Fit("AGaus","Q0RW");
         TF1 *Novosib = new TF1("Novosib", this, &MpdPidQA::Novosibirsk, XFUNCMIN, XFUNCMAX, 4, "MpdPidQA",
   "Novosibirsk"); Novosib->SetParameters(Gaus->GetParameter(0), 0.01, Gaus->GetParameter(2), Gaus->GetParameter(1));
         it->second.dEdXPart[i]->Fit("Novosib","Q0RW");

         intF_MultX = 0.0; intF = 0.0;
         for ( TF1* BBF : vecParBB ) {
            /// Xlow and Xhigh are within one Bethe-Bloch function
            if ( ( BBF->GetXmin() <= Xlow[i+ibegloc] ) && ( BBF->GetXmax() > Xlow[i+ibegloc] ) && ( BBF->GetXmin() <=
   Xhigh[i+ibegloc] ) && ( BBF->GetXmax() > Xhigh[i+ibegloc] ) ) { for ( TF1* BBF_MultX : vecParBBMultX ) { FName =
   BBF->GetName(); FName.Append("_MultX"); FMXName = BBF_MultX->GetName(); if ( FMXName == FName ) { intF_MultX =
   BBF_MultX->Integral( Xlow[i+ibegloc], Xhigh[i+ibegloc] ); intF = BBF->Integral( Xlow[i+ibegloc], Xhigh[i+ibegloc] );
                     mom = intF_MultX / intF;
                  }
               }
            /// Xlow and Xhigh are within two functions
            } else if ( ( BBF->GetXmin() <= Xlow[i+ibegloc] ) && ( BBF->GetXmax() > Xlow[i+ibegloc] ) ) {
               intF += BBF->Integral( Xlow[i+ibegloc], BBF->GetXmax() );
               for ( TF1* BBF_MultX : vecParBBMultX ) {
                  FName = BBF->GetName(); FName.Append("_MultX");
                  FMXName = BBF_MultX->GetName();
                  if ( FMXName == FName ) {
                     intF_MultX += BBF_MultX->Integral( Xlow[i+ibegloc], BBF_MultX->GetXmax() );
                  }
               }
            } else if ( ( BBF->GetXmin() > Xlow[i+ibegloc] ) && ( BBF->GetXmax() < Xhigh[i+ibegloc] ) ) {
               if ( intF != 0.0 ) {
                  intF += BBF->Integral( BBF->GetXmin(), BBF->GetXmax() );
                  for ( TF1* BBF_MultX : vecParBBMultX ) {
                     FName = BBF->GetName(); FName.Append("_MultX");
                     FMXName = BBF_MultX->GetName();
                     if ( FMXName == FName ) {
                        intF_MultX += BBF_MultX->Integral( BBF->GetXmin(), BBF->GetXmax() );
                     }
                  }
               }
            } else if ( ( BBF->GetXmin() <= Xhigh[i+ibegloc] ) && ( BBF->GetXmax() > Xhigh[i+ibegloc] ) ) {
               if ( intF != 0.0 ) {
                  intF += BBF->Integral( BBF->GetXmin(), Xhigh[i+ibegloc] );
                  for ( TF1* BBF_MultX : vecParBBMultX ) {
                     FName = BBF->GetName(); FName.Append("_MultX");
                     FMXName = BBF_MultX->GetName();
                     if ( FMXName == FName ) {
                        intF_MultX += BBF_MultX->Integral( BBF->GetXmin(), Xhigh[i+ibegloc] );
                        mom = intF_MultX / intF;
                     }
                  }
               }
            } else mom = -999.0;
         }

         Ygraph[i] = ( fTrackingState == MpdPidUtils::kCFHM ) ? AGaus->GetParameter(1) : Novosib->GetParameter(3);
         Ygraph[i] /= GetDedxParam(mom,it->first);
         sigma1[i] = AGaus->GetParameter(2); sigma2[i] = (1. + AGaus->GetParameter(3)) * AGaus->GetParameter(2);
         sigma1[i] /= GetDedxParam(mom,it->first); sigma2[i] /= GetDedxParam(mom,it->first);
         delete Gaus; delete AGaus; delete Novosib;
      }
      TGraphAsymmErrors *gr = new TGraphAsymmErrors(nQAHistsloc,Xgraph,Ygraph,Xerrgraph,Xerrgraph,sigma1,sigma2);
      FName = "Graph_"; FName.Append( MpdPidUtils::cParticleShortName[it->first] );
      gr->SetName(FName);
      graphs.insert( pair<MpdPidUtils::ePartType,TGraphAsymmErrors*>(it->first,gr) );
      delete Xgraph; delete Ygraph; delete Xerrgraph; delete sigma;
   }
*/
   if (!dir.EndsWith("/")) dir.Append("/");
   TFile *fFile = new TFile(dir + "dEdxHists.root", "RECREATE");

   for (auto it = fEnLossMap.begin(); it != fEnLossMap.end(); ++it) {
      for (Int_t i = 0; i < (it->second.iend - it->second.ibeg + 1); i++) {
         it->second.dEdXPart[i]->Write();
      }
      it->second.BetheBlochHist->Write();
   }
   // for ( auto it = graphs.begin(); it != graphs.end(); ++it ) {
   //	it->second->Write();
   // }
   fSumBetheBlochHists->Write();
   fChBetheBlochHists->Write();
   fFile->Close();
}

void MpdPidQA::Getm2QA(TString dir)
{
   if (!dir.EndsWith("/")) dir.Append("/");
   TFile *fFile = new TFile(dir + "m2Hists.root", "RECREATE");

   for (auto it = fMSquaredMap.begin(); it != fMSquaredMap.end(); ++it) {
      it->second->Write();
   }
   fm2LightHist->Write();
   fm2HeavyHist->Write();
   fFile->Close();
}

void MpdPidQA::GetAmplQA(TString dir)
{
   if (!dir.EndsWith("/")) dir.Append("/");
   TFile *fFile = new TFile(dir + "ParticleYields.root", "RECREATE");

   for (auto it = fAbundanceMap.begin(); it != fAbundanceMap.end(); ++it) {
      for (TH1D *hYield : it->second) {
         hYield->Write();
      }
   }
   fFile->Close();
}

void MpdPidQA::GetEffContQA(TString dir)
{
   if (!dir.EndsWith("/")) dir.Append("/");
   TFile *fFile = new TFile(dir + "EffCont.root", "RECREATE");

   for (auto it = fEffContMap.begin(); it != fEffContMap.end(); ++it) {
      for (vecTH1Dptrs vecEC : it->second) {
         vecEC[2]->Divide(vecEC[0]);
         vecEC[3]->Divide(vecEC[1]);
         vecEC[2]->Write();
         vecEC[3]->Write();
      }
   }
   fFile->Close();
}

void MpdPidQA::GetEffContQA_Raw(TString dir)
{
   if (!dir.EndsWith("/")) dir.Append("/");
   TFile *fFile = new TFile(dir + "EffCont.root", "RECREATE");

   for (auto it = fEffContMap.begin(); it != fEffContMap.end(); ++it) {
      for (vecTH1Dptrs vecEC : it->second) {
         vecEC[0]->Write();
         vecEC[1]->Write();
         vecEC[2]->Write();
         vecEC[3]->Write();
      }
   }
   fFile->Close();
}

Double_t MpdPidQA::Novosibirsk(Double_t *x, Double_t *par)
{
   Double_t tail  = par[1];
   Double_t width = par[2];
   Double_t peak  = par[3];
   if (TMath::Abs(tail) < 1.e-7) return par[0] * TMath::Exp(-0.5 * TMath::Power(((x[0] - peak) / width), 2));
   Double_t arg = 1.0 - (x[0] - peak) * tail / width;
   if (arg < 1.e-6) return 0.0; /// Argument of logaritem negative. Real continuation -> function equals zero

   Double_t              log         = TMath::Log(arg);
   static const Double_t xi          = 2.3548200450309494; /// 2 Sqrt( Ln(4) )
   Double_t              width_zero  = (2.0 / xi) * TMath::ASinH(tail * xi * 0.5);
   Double_t              width_zero2 = width_zero * width_zero;
   Double_t              exponent    = (-0.5 / (width_zero2)*log * log) - (width_zero2 * 0.5);

   return par[0] * TMath::Exp(exponent);
}

ClassImp(MpdPidQA);
