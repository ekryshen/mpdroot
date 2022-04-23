#include <iostream>
#include <fstream> // std::ifstream

#include "MpdKalmanFilter.h"
#include "MpdMCTrack.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdEmcClusterKI.h"
#include "MpdParticle.h"
#include "MpdVertex.h"
#include "MpdEvent.h"
#include "MpdEmcGeoUtils.h"

#include "MpdConvPi0.h"
#include "TFile.h"

ClassImp(MpdConvPi0);

MpdConvPi0::MpdConvPi0(const char *name, const char *outputName) : MpdAnalysisTask(name, outputName)
{
   mParamConfig = outputName;
}

void MpdConvPi0::UserInit()
{

   mParams.ReadFromFile(mParamConfig);
   mParams.Print();

   // // Setup branches
   // mChain->SetBranchAddress("MCEventHeader.", &mMCHeader);
   // mChain->SetBranchAddress("MCTrack", &mMCTracks);
   // mChain->SetBranchAddress("MPDEvent.", &mMPDEvent);

   // mEMCClusters  = new TObjArray() ;
   // mChain->SetBranchAddress("EmcCluster",&mEMCClusters);
   // mChain->SetBranchAddress("Vertex",&mVertexes);

   // mChain->SetBranchAddress("TpcKalmanTrack",&mKalmanTracks);

   // Prepare histograms etc.
   fOutputList = new TList();
   fOutputList->SetOwner(kTRUE);

   TH1::AddDirectory(kFALSE); // sets a global switch disabling the reference to histos in gROOT and their overwriting

   // General QA
   mhEvents = new TH1F("hEvents", "Number of events", 10, 0., 10.);
   fOutputList->Add(mhEvents);
   mhVertex = new TH1F("hVertex", "Event vertex distribution", 50, -100., 100.);
   fOutputList->Add(mhVertex);
   mhCentrality = new TH1F("hCentrality", "Centrality distribution", 100, 0., 100.);
   fOutputList->Add(mhCentrality);

   // V0selection
   const int   nPtbin = 100;
   const float pTmin  = 0.;
   const float pTmax  = 5.;
   mhCutEff           = new TH2F("hCutEff", "track cut efficiency", 21, 0., 21., nPtbin, pTmin, pTmax);
   fOutputList->Add(mhCutEff);
   mhNhits = new TH1F("hNhits", "Number of hits per track", 100, 0., 100.);
   fOutputList->Add(mhNhits);
   mhTracks = new TH2F("hTracks", "track occupancy pt vs eta", nPtbin, pTmin, pTmax, 100, -1.5, 1.5);
   fOutputList->Add(mhTracks);
   mhProbEl = new TH2F("hProbEl", "Electron probability", 100, 0., 1., nPtbin, pTmin, pTmax);
   fOutputList->Add(mhProbEl);
   mhdEdx = new TH2F("hdEdx", "dEdx", 100, -10, 10., nPtbin, pTmin, pTmax);
   fOutputList->Add(mhdEdx);
   mhConvMap = new TH3F("hConvMap", "Conversion map (r,phi,z)", 100, 0., 280., 100, 0., TMath::Pi(), 100, -200., 200.);
   fOutputList->Add(mhConvMap);
   mhAlpha = new TH2F("hAlpha", "#alpha distribution", 100, 0., 1., nPtbin, pTmin, pTmax);
   fOutputList->Add(mhAlpha);
   mhChi2 = new TH2F("hChi2", "#chi^{2}", 100, 0., 50., nPtbin, pTmin, pTmax);
   fOutputList->Add(mhChi2);
   mhDist = new TH2F("hDist", "track DCA", 100, 0., 10., nPtbin, pTmin, pTmax);
   fOutputList->Add(mhDist);
   hmMassEE = new TH2F("mEE", "m_{ee}", 100, 0., 0.3, nPtbin, pTmin, pTmax);
   fOutputList->Add(hmMassEE);
   mhCosPsi = new TH2F("cosPsi", "cos(#psi)", 100, 0., TMath::Pi(), nPtbin, pTmin, pTmax);
   fOutputList->Add(mhCosPsi);
   mhArmPo = new TH2F("Armenteros", "Armenteros", 100, -1, 1, 100, 0, 0.3);
   fOutputList->Add(mhArmPo);
   mhConvSp = new TH2F("ConvSp", "Conv Sp", nPtbin, pTmin, pTmax, 100, -1.5, 1.5);
   fOutputList->Add(mhConvSp);
   mhAsym = new TH2F("Asymetry", "Asymetry", 200, 0, 1, nPtbin, pTmin, pTmax);
   fOutputList->Add(mhAsym);
   if (isMC) {
      mhProbElTrue = new TH2F("hProbElTrue", "Electron probability, true electrons", 100, 0., 1., nPtbin, pTmin, pTmax);
      fOutputList->Add(mhProbElTrue);
      mhdEdxTrue = new TH2F("hdEdxTrue", "dEdx", 100, -10, 10., nPtbin, pTmin, pTmax);
      fOutputList->Add(mhdEdxTrue);
      mhConvMapTrue =
         new TH3F("hConvMapTrue", "Conversion map (r,phi,z)", 100, 0., 280., 100, 0., TMath::Pi(), 100, -200., 200.);
      fOutputList->Add(mhConvMapTrue);
      mhAlphaTrue = new TH2F("hAlphaTrue", "#alpha distribution", 100, 0., 1., nPtbin, pTmin, pTmax);
      fOutputList->Add(mhAlphaTrue);
      mhChi2True = new TH2F("hChi2True", "#chi^{2}", 100, 0., 50., nPtbin, pTmin, pTmax);
      fOutputList->Add(mhChi2True);
      mhDistTrue = new TH2F("hDistTrue", "track DCA", 100, 0., 10., nPtbin, pTmin, pTmax);
      fOutputList->Add(mhDistTrue);
      hmMassEETrue = new TH2F("mEETrue", "m_{ee}", 100, 0., 0.3, nPtbin, pTmin, pTmax);
      fOutputList->Add(hmMassEETrue);
      mhCosPsiTrue = new TH2F("cosPsiTrue", "cos(#psi)", 100, 0., TMath::Pi(), nPtbin, pTmin, pTmax);
      fOutputList->Add(mhCosPsiTrue);
      mhArmPoTrue = new TH2F("ArmenterosTrue", "Armenteros", 100, -1, 1, 100, 0, 0.3);
      fOutputList->Add(mhArmPoTrue);
      mhAsymTrue = new TH2F("AsymetryTrue", "Asymetry", 200, 0, 1, nPtbin, pTmin, pTmax);
      fOutputList->Add(mhAsymTrue);
      mhConvSpTrue = new TH2F("ConvSpTrue", "Conv Sp", nPtbin, pTmin, pTmax, 100, -1.5, 1.5);
      fOutputList->Add(mhConvSpTrue);
   }

   // Cluster selection
   mhCluCutEff = new TH2F("hCluCutEff", "cluster cut efficiency", 10, 0., 10., nPtbin, pTmin, pTmax);
   fOutputList->Add(mhCluCutEff);
   // Inv mass histos
   const int   nMbins = 200;
   const float mMax   = 1.;
   for (int cen = 0; cen < mHistoCentBins; cen++) {
      mhRealCalo[cen] = new TH2F(Form("hRealCalo_cen%d", cen), Form("Real inv mass, calorimeter"), nMbins, 0., mMax,
                                 nPtbin, pTmin, pTmax);
      fOutputList->Add(mhRealCalo[cen]);
      mhMixedCalo[cen] = new TH2F(Form("hMixedCalo_cen%d", cen), Form("Mixed inv mass, calorimeter"), nMbins, 0., mMax,
                                  nPtbin, pTmin, pTmax);
      fOutputList->Add(mhMixedCalo[cen]);
      mhRealHybrid[cen] = new TH2F(Form("hRealHybrid_cen%d", cen), Form("Real inv mass, hybrid"), nMbins, 0., mMax,
                                   nPtbin, pTmin, pTmax);
      fOutputList->Add(mhRealHybrid[cen]);
      mhMixedHybrid[cen] = new TH2F(Form("hMixedHybrid_cen%d", cen), Form("Mixed inv mass, hybrid"), nMbins, 0., mMax,
                                    nPtbin, pTmin, pTmax);
      fOutputList->Add(mhMixedHybrid[cen]);
      mhRealConv[cen] = new TH2F(Form("hRealConv_cen%d", cen), Form("Real inv mass, conversion"), nMbins, 0., mMax,
                                 nPtbin, pTmin, pTmax);
      fOutputList->Add(mhRealConv[cen]);
      mhMixedConv[cen] = new TH2F(Form("hMixedConversion_cen%d", cen), Form("Mixed inv mass, conversion"), nMbins, 0.,
                                  mMax, nPtbin, pTmin, pTmax);
      fOutputList->Add(mhMixedConv[cen]);
      if (isMC) {
         mhRealCaloTrue[cen] = new TH2F(Form("hRealCaloTrue_cen%d", cen), Form("Real inv mass, calorimeter"), nMbins,
                                        0., mMax, nPtbin, pTmin, pTmax);
         fOutputList->Add(mhRealCaloTrue[cen]);
         mhRealHybridTrue[cen] = new TH2F(Form("hRealHybridTrue_cen%d", cen), Form("Real inv mass, hybrid"), nMbins, 0.,
                                          mMax, nPtbin, pTmin, pTmax);
         fOutputList->Add(mhRealHybridTrue[cen]);
         mhRealConvTrue[cen] = new TH2F(Form("hRealConvTrue_cen%d", cen), Form("Real inv mass, conversion"), nMbins, 0.,
                                        mMax, nPtbin, pTmin, pTmax);
         fOutputList->Add(mhRealConvTrue[cen]);

         mhRealCaloTrueEl[cen] = new TH2F(Form("hRealCaloTrueEl_cen%d", cen), Form("Real inv mass, calorimeter"),
                                          nMbins, 0., mMax, nPtbin, pTmin, pTmax);
         fOutputList->Add(mhRealCaloTrueEl[cen]);
         mhRealHybridTrueEl[cen] = new TH2F(Form("hRealHybridTrueEl_cen%d", cen), Form("Real inv mass, hybrid"), nMbins,
                                            0., mMax, nPtbin, pTmin, pTmax);
         fOutputList->Add(mhRealHybridTrueEl[cen]);
         mhRealConvTrueEl[cen] = new TH2F(Form("hRealConvTrueEl_cen%d", cen), Form("Real inv mass, conversion"), nMbins,
                                          0., mMax, nPtbin, pTmin, pTmax);
         fOutputList->Add(mhRealConvTrueEl[cen]);

         mhRealCaloTrueAll[cen] = new TH2F(Form("hRealCaloTrueAll_cen%d", cen), Form("Real inv mass, calorimeter"),
                                           nMbins, 0., mMax, nPtbin, pTmin, pTmax);
         fOutputList->Add(mhRealCaloTrueAll[cen]);
         mhRealHybridTrueAll[cen] = new TH2F(Form("hRealHybridTrueAll_cen%d", cen), Form("Real inv mass, hybrid"),
                                             nMbins, 0., mMax, nPtbin, pTmin, pTmax);
         fOutputList->Add(mhRealHybridTrueAll[cen]);
         mhRealConvTrueAll[cen] = new TH2F(Form("hRealConvTrueAll_cen%d", cen), Form("Real inv mass, conversion"),
                                           nMbins, 0., mMax, nPtbin, pTmin, pTmax);
         fOutputList->Add(mhRealConvTrueAll[cen]);
      }
   }
   if (isMC) {
      // spectra of primary pi0.eta,...
      for (int cen = 0; cen < mHistoCentBins; cen++) {
         hPrimPi0[cen] = new TH2F(Form("PrimaryPi0_cen%d", cen), "Centrality", nPtbin, pTmin, pTmax, 100, -2., 2.);
         fOutputList->Add(hPrimPi0[cen]);
      }
   }
}
//--------------------------------------
void MpdConvPi0::ProcessEvent(MpdAnalysisEvent &event)
{

   if (!isInitialized) {
      mKF = MpdKalmanFilter::Instance();
      mKHit.SetType(MpdKalmanHit::kFixedR);
      mPID          = new MpdPid(mParams.mPIDsigM, mParams.mPIDsigE, mParams.mPIDenergy, mParams.mPIDkoeff,
                                 mParams.mPIDgenerator, mParams.mPIDtracking, mParams.mPIDparticles);
      isInitialized = true;
   }

   if (!selectEvent(event)) { ///
      return;
   }

   mEMCClusters  = event.fEMCCluster;
   mKalmanTracks = event.fTPCKalmanTrack;
   if (isMC) {
      mMCTracks = event.fMCTrack;
      for (int i = 0; i < mMCTracks->GetEntriesFast(); i++) {
         MpdMCTrack *pr = (static_cast<MpdMCTrack *>(mMCTracks->At(i)));
         if (pr->GetPdgCode() == 111) {
            if (pr->GetStartX() * pr->GetStartX() + pr->GetStartY() * pr->GetStartY() < 1.) {
               TVector3 momentum;
               pr->GetMomentum(momentum);
               hPrimPi0[mCenBin]->Fill(momentum.Pt(), momentum.Y());
            }
         }
      }
   }

   selectConversion(event);

   selectClusters(event);

   processHistograms(event);
}

void MpdConvPi0::Finish()
{
   // Post-scan processing not needed
}

//--------------------------------------
bool MpdConvPi0::selectEvent(MpdAnalysisEvent &event)
{

   mhEvents->Fill(0.5);
   // first test if event filled?
   if (!event.fVertex) { // if even vertex not filled, skip event
      return false;
   }
   // Vertex z coordinate
   MpdVertex *vertex = (MpdVertex *)event.fVertex->First();
   vertex->Position(mPrimaryVertex);
   mhVertex->Fill(mPrimaryVertex.Z());
   if (fabs(mPrimaryVertex.Z()) > mParams.mZvtxCut) {
      return false;
   }
   mZvtxBin = 0.5 * (mPrimaryVertex.Z() / mParams.mZvtxCut + 1) * nMixEventZ;
   if (mZvtxBin < 0) mZvtxBin = 0;
   if (mZvtxBin >= nMixEventZ) mZvtxBin = nMixEventZ - 1;
   mhEvents->Fill(1.5);

   // Centrality
   //  mhCentrality->Fill(cen) ;
   mCenBin = 0;
   mhEvents->Fill(2.5);

   // ZCD vs TPC (pileup?)
   mhEvents->Fill(3.5);

   // Eventplane  TODO
   mRPBin = 0;
   mhEvents->Fill(4.5);

   return true;
}
//--------------------------------------
void MpdConvPi0::selectConversion(MpdAnalysisEvent &event)
{

   // Select V0 in this event
   mV0.clear();

   mMpdGlobalTracks = event.fMPDEvent->GetGlobalTracks();
   int       ntr    = mMpdGlobalTracks->GetEntriesFast();
   MpdPhoton v;
   for (long int i = 0; i < ntr - 1; i++) {
      MpdTrack          *mpdtrack1 = (MpdTrack *)mMpdGlobalTracks->UncheckedAt(i);
      MpdTpcKalmanTrack *tr1       = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(i);
      if (!selectTrack(mpdtrack1)) {
         continue;
      }
      for (long int j = i + 1; j < ntr; j++) {
         MpdTrack          *mpdtrack2 = (MpdTrack *)mMpdGlobalTracks->UncheckedAt(j);
         MpdTpcKalmanTrack *tr2       = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(j);
         if (!selectTrack(mpdtrack2)) {
            continue;
         }
         if (createSelectV0(mpdtrack1, tr1, mpdtrack2, tr2, v)) {
            mV0.emplace_back(v.X(), v.Y(), v.Z(), v.T());
            mV0.back().setPrimary(v.primary());
            mV0.back().setTr1(i);
            mV0.back().setTr2(j);
         }
      }
   }
}
//--------------------------------------
void MpdConvPi0::selectClusters(MpdAnalysisEvent &event)
{
   // Select EMC clusters in event
   const Float_t par0[2] = {-7.02851e-005, -7.76476e-002};
   const Float_t par1[2] = {-4.45930e-005, -1.01164e-002};

   mClusters.clear();
   int n = mEMCClusters->GetEntriesFast();
   for (int i = n; i--;) {
      MpdEmcClusterKI *clu = (MpdEmcClusterKI *)mEMCClusters->At(i);

      float e = Nonlinearity(clu->GetE());
      mhCluCutEff->Fill(0., e);
      if (e < mParams.mCluEmin) {
         continue;
      }

      if (clu->GetMultiplicity() < mParams.mCluMult) {
         continue;
      }
      mhCluCutEff->Fill(1., e);

      double dx   = clu->GetX() - mPrimaryVertex.X();
      double dy   = clu->GetY() - mPrimaryVertex.Y();
      double dz   = clu->GetZ() - mPrimaryVertex.Z();
      double r    = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
      double time = clu->GetTime() - r / 29979245800. * 1.e+9; // Time in ns
      if (fabs(tofCut(time, e)) > mParams.mCluTof) {
         continue;
      }
      mhCluCutEff->Fill(2., e);

      float l1, l2; // Dispersion axis
      clu->GetLambdas(l1, l2);
      // correct for z
      float izMax = TMath::Abs(32. * clu->GetZ());
      l1          = l1 * 1.6 / (1.6 + 0.0002 * izMax * izMax);
      l2          = l2 * 3.2 / (3.2 + 0.000023 * izMax * izMax * izMax);

      //    if(e>mParams.mCluDispEmin &&  lambdaCut(l1,l2,e) > mParams.mCluDisp){
      //      continue ;
      //    }
      mhCluCutEff->Fill(3., e);

      float pp0      = par0[0] + par0[1] * mPrimaryVertex.Z();
      float pp1      = par1[0] + par1[1] * mPrimaryVertex.Z();
      float z_shift1 = pp0 + pp1 * log(e);
      if (distCPV(clu->GetDPhi(), clu->GetDZ() + z_shift1, e) < mParams.mCluCPV) {
         continue;
      }
      mhCluCutEff->Fill(4., e);

      mClusters.emplace_back(dx / r * e, dy / r * e, dz / r * e, e);
      mClusters.back().setTr1(i);
      if (clu->GetNumberOfTracks() > 0) {
         int   trackId;
         float edep;
         clu->GetMCTrack(0, trackId, edep);
         mClusters.back().setPrimary(trackId);
      }
   }
}

void MpdConvPi0::processHistograms(MpdAnalysisEvent &event)
{
   // Fill Real, Mixed distributions and update mixed array

   // Real
   int nClu = mClusters.size();
   int nV0  = mV0.size();
   for (int i = 0; i < nClu - 1; i++) {
      for (int j = i + 1; j < nClu; j++) {
         TLorentzVector sum = mClusters[i] + mClusters[j];
         mhRealCalo[mCenBin]->Fill(sum.M(), sum.Pt());
         long int ip = IsSameParent(mClusters[i].primary(), mClusters[j].primary());
         if (ip >= 0) {
            mhRealCaloTrueAll[mCenBin]->Fill(sum.M(), sum.Pt());
            if (static_cast<MpdMCTrack *>(mMCTracks->At(ip))->GetPdgCode() == 111) {
               mhRealCaloTrue[mCenBin]->Fill(sum.M(), sum.Pt());
            }
            if (TMath::Abs(static_cast<MpdMCTrack *>(mMCTracks->At(ip))->GetPdgCode()) == 11 ||
                static_cast<MpdMCTrack *>(mMCTracks->At(ip))->GetPdgCode() == 22) {
               mhRealCaloTrueEl[mCenBin]->Fill(sum.M(), sum.Pt());
            }
         }
      }
   }

   for (int i = 0; i < nClu; i++) {
      for (int j = 0; j < nV0; j++) {
         if (TestHybrid(mClusters[i], mV0[j])) {
            TLorentzVector sum = mClusters[i] + mV0[j];
            mhRealHybrid[mCenBin]->Fill(sum.M(), sum.Pt());
            long int ip = IsSameParent(mClusters[i].primary(), mV0[j].primary());
            if (ip >= 0) {
               mhRealHybridTrueAll[mCenBin]->Fill(sum.M(), sum.Pt());
               if (static_cast<MpdMCTrack *>(mMCTracks->At(ip))->GetPdgCode() == 111) {
                  mhRealHybridTrue[mCenBin]->Fill(sum.M(), sum.Pt());
               }
               if (TMath::Abs(static_cast<MpdMCTrack *>(mMCTracks->At(ip))->GetPdgCode()) == 11 ||
                   static_cast<MpdMCTrack *>(mMCTracks->At(ip))->GetPdgCode() == 22) {
                  mhRealHybridTrueEl[mCenBin]->Fill(sum.M(), sum.Pt());
               }
            }
         }
      }
   }
   for (int i = 0; i < nV0 - 1; i++) {
      for (int j = i + 1; j < nV0; j++) {
         TLorentzVector sum = mV0[i] + mV0[j];
         mhRealConv[mCenBin]->Fill(sum.M(), sum.Pt());
         long int ip = IsSameParent(mV0[i].primary(), mV0[j].primary());
         if (ip >= 0) {
            mhRealConvTrueAll[mCenBin]->Fill(sum.M(), sum.Pt());
            if (static_cast<MpdMCTrack *>(mMCTracks->At(ip))->GetPdgCode() == 111) {
               mhRealConvTrue[mCenBin]->Fill(sum.M(), sum.Pt());
            }
            if (TMath::Abs(static_cast<MpdMCTrack *>(mMCTracks->At(ip))->GetPdgCode()) == 11 ||
                static_cast<MpdMCTrack *>(mMCTracks->At(ip))->GetPdgCode() == 22) {
               mhRealConvTrueEl[mCenBin]->Fill(sum.M(), sum.Pt());
            }
         }
      }
   }

   // Mixed
   // calculate bin from zVertex-centrality-reaction plane
   int mixBin = mZvtxBin * nMixEventCent * nMixEventRP + mCenBin * nMixEventRP + mRPBin;

   for (auto &vm : mMixClu[mixBin]) {
      for (auto &v : mClusters) {
         TLorentzVector sum = v + vm;
         mhMixedCalo[mCenBin]->Fill(sum.M(), sum.Pt());
      }
   }
   for (auto &vm : mMixV0[mixBin]) {
      for (auto &v : mClusters) {
         TLorentzVector sum = v + vm;
         mhMixedHybrid[mCenBin]->Fill(sum.M(), sum.Pt());
      }
   }
   for (auto &vm : mMixClu[mixBin]) {
      for (auto &v : mV0) {
         TLorentzVector sum = v + vm;
         mhMixedHybrid[mCenBin]->Fill(sum.M(), sum.Pt());
      }
   }
   for (auto &vm : mMixV0[mixBin]) {
      for (auto &v : mV0) {
         TLorentzVector sum = v + vm;
         mhMixedConv[mCenBin]->Fill(sum.M(), sum.Pt());
      }
   }

   // Append new particles to queue and remove those at the beginning
   for (auto &v : mV0) {
      mMixV0[mixBin].emplace_back(v);
   }
   while (mMixV0[mixBin].size() > (UInt_t)mMaxMixSize) {
      mMixV0[mixBin].pop_front();
   }
   for (auto &v : mClusters) {
      mMixClu[mCenBin].emplace_back(v);
   }
   while (mMixClu[mixBin].size() > (UInt_t)mMaxMixSize) {
      mMixClu[mixBin].pop_front();
   }
}

bool MpdConvPi0::selectTrack(MpdTrack *mpdtrack)
{

   mhNhits->Fill(mpdtrack->GetNofHits());
   float pt = TMath::Abs(mpdtrack->GetPt());
   mhCutEff->Fill(0., pt);
   if (mpdtrack->GetNofHits() < mParams.mNofHitsCut) return false; // nhits > 10
   mhCutEff->Fill(1., pt);
   mhTracks->Fill(pt, mpdtrack->GetEta());
   if (fabs(mpdtrack->GetEta()) > mParams.mEtaCut) return false; //|eta| < 1.0
   if (pt < mParams.mPtminCut) return false;                     // pT > 50 MeV/c
   mhCutEff->Fill(2., pt);

   int charge;
   if (mpdtrack->GetPt() < 0)
      charge = 1;
   else
      charge = -1;

   bool isGoodPID;
   if (mpdtrack->GetTofFlag() == 2 || mpdtrack->GetTofFlag() == 6) {
      isGoodPID =
         mPID->FillProbs(pt * TMath::CosH(mpdtrack->GetEta()), mpdtrack->GetdEdXTPC(), mpdtrack->GetTofMass2(), charge);
   } else {
      isGoodPID = mPID->FillProbs(pt * TMath::CosH(mpdtrack->GetEta()), mpdtrack->GetdEdXTPC(), charge);
   }

   mhProbEl->Fill(mPID->GetProbEl(), fabs(mpdtrack->GetPt()));
   bool isElectron = false;
   if (isMC) { // same for true electron tracks
      long int prim1 = mpdtrack->GetID();
      if (prim1 >= 0) {
         isElectron = abs((static_cast<MpdMCTrack *>(mMCTracks->At(prim1)))->GetPdgCode()) == 11;
      }
   }
   if (isElectron) {
      mhProbElTrue->Fill(mPID->GetProbEl(), pt);
   }

   if (isGoodPID && mPID->GetProbEl() < mParams.mProbElCut) {
      return false;
   }
   mhCutEff->Fill(3., pt);
   float dEdx = dEdx_sigma(mpdtrack->GetdEdXTPC(), sqrt(pow(pt, 2) + pow(mpdtrack->GetPz(), 2)));
   mhdEdx->Fill(dEdx, pt);
   if (isElectron) { // same for true electrontracks
      mhdEdxTrue->Fill(dEdx, pt);
   }
   if ((fabs(dEdx) < mParams.mdEdxSigmaCut) &&
       (fabs(Beta_sigma(mpdtrack->GetTofBeta(), sqrt(pow(pt, 2) + pow(mpdtrack->GetPz(), 2)))) <
           mParams.mBetaSigmaCut &&
        (mpdtrack->GetTofFlag() == 2 || mpdtrack->GetTofFlag() == 6))) {
      mhCutEff->Fill(4., pt);
      return true;
   } else {
      return false;
   }
   return false;
}

bool MpdConvPi0::createSelectV0(MpdTrack *tr1, MpdTpcKalmanTrack *ktr1, MpdTrack *tr2, MpdTpcKalmanTrack *ktr2,
                                MpdPhoton &v)
{
   // Construct and check V0

   // Use opposite charge tracks
   int charge1, charge2;
   if (tr1->GetPt() < 0)
      charge1 = 1;
   else
      charge1 = -1;
   if (tr2->GetPt() < 0)
      charge2 = 1;
   else
      charge2 = -1;

   // reject same sign pairs
   if (charge1 * charge2 > 0) return false;

   // Will be used for extrapolation
   MpdTpcKalmanTrack trCorK1(*ktr1);
   MpdHelix          helix1 = MakeHelix(trCorK1);
   MpdParticle       el1(trCorK1, 0);
   if (charge1 > 0) el1.SetPdg(-11);
   if (charge1 < 0) el1.SetPdg(11);
   el1.SetMass();

   MpdTpcKalmanTrack trCorK2(*ktr2);
   MpdHelix          helix2 = MakeHelix(trCorK2);
   MpdParticle       el2(trCorK2, 0);

   if (charge2 > 0) el2.SetPdg(-11);
   if (charge2 < 0) el2.SetPdg(11);
   el2.SetMass();

   // pair
   mPartK.clear();
   mPartK.emplace_back(&el1);
   mPartK.emplace_back(&el2);

   MpdParticle gamEE;
   float       chi2 = TMath::Abs(gamEE.BuildMother(mPartK));
   float       pt   = gamEE.Pt();
   mhChi2->Fill(chi2, pt);
   mhCutEff->Fill(10., pt);

   if (pt < 0.005) { // to avoid fpe
      return false;
   }
   mhCutEff->Fill(11., pt);

   bool     isTrue         = false; // is true conv pair?
   long int commonParentId = -1;
   if (isMC) { // same for true electrontracks
      long int prim1 = tr1->GetID();
      long int prim2 = tr2->GetID();
      commonParentId = IsSameParent(prim1, prim2);
      if (commonParentId >= 0) { // there is common parent
         isTrue = (static_cast<MpdMCTrack *>(mMCTracks->At(commonParentId))->GetPdgCode() == 22);
      }
   }
   if (isTrue) {
      mhChi2True->Fill(chi2, pt);
   }
   if (chi2 > mParams.mChi2Cut) {
      return false;
   }
   mhCutEff->Fill(12., pt);

   TVector3 v0(gamEE.Getx()(0, 0), gamEE.Getx()(1, 0), gamEE.Getx()(2, 0));
   v0 -= mPrimaryVertex;
   // float decay = v0.Mag();

   mhConvMap->Fill(gamEE.Getx()(0, 0), gamEE.Getx()(1, 0), gamEE.Getx()(2, 0));
   if (isTrue) {
      mhConvMapTrue->Fill(gamEE.Getx()(0, 0), gamEE.Getx()(1, 0), gamEE.Getx()(2, 0));
   }

   float rConv = TMath::Sqrt(pow(gamEE.Getx()(0, 0), 2) + pow(gamEE.Getx()(1, 0), 2));
   if (rConv < mParams.mMinR2Cut || rConv > mParams.mMaxR2Cut) {
      return false;
   }
   mhCutEff->Fill(13., pt);

   float angle; // disth,
   angle = v0.Angle(gamEE.Momentum3());
   mhAlpha->Fill(angle, pt);
   if (isTrue) { // same for true electrontracks
      mhAlphaTrue->Fill(angle, pt);
   }
   if (angle > mParams.mAlphaCut) {
      return false;
   }

   mhCutEff->Fill(14., pt);

   // if( ePos->R() <= ((TMath::Abs(ePos->Vz()) * fLineCutZRSlope) - fLineCutZValue)){
   //   return kFALSE;  // line cut to exclude regions where we do not reconstruct
   // } else if ( fEtaCutMin != -0.1 &&   ePos->R() >= ((TMath::Abs(ePos->Vz()) * fLineCutZRSlopeMin) -
   // fLineCutZValueMin)){
   //   return kFALSE;
   // }

   std::pair<float, float> paths = helix1.pathLengths(helix2);
   TVector3                p1    = helix1.at(paths.first);
   TVector3                p2    = helix2.at(paths.second);
   p1 -= p2;
   float dist = p1.Mag(); // Closest distance between daughters

   mhDist->Fill(dist, pt);
   if (isTrue) { // same for true electrontracks
      mhDistTrue->Fill(dist, pt);
   }
   if (dist > mParams.mDistCut) {
      return false;
   }
   mhCutEff->Fill(15., pt);

   hmMassEE->Fill(gamEE.GetMass(), pt);
   if (isTrue) { // same for true electrontracks
      hmMassEETrue->Fill(gamEE.GetMass(), pt);
   }
   if (gamEE.GetMass() > mParams.mMassCut) {
      return false;
   }
   mhCutEff->Fill(16., pt);

   // Pair_chi_1[n_ks] = el1.Chi2Vertex(vertex);
   // Pair_chi_2[n_ks] = el2.Chi2Vertex(vertex);
   // A-P cut
   //  Gamma selection based on QT from Armenteros
   // propagate trCorK1,trCorK2 to conversion point

   MpdKalmanHit hitTmp;
   hitTmp.SetType(MpdKalmanHit::kFixedR);
   hitTmp.SetPos(trCorK1.GetPos());
   trCorK1.SetParamNew(*trCorK1.GetParam());
   trCorK1.SetPos(trCorK1.GetPosNew());
   trCorK1.ReSetWeight();
   //  TMatrixDSym w = *trCorK1.GetWeight(); // save current weight matrix
   mKHit.SetPos(rConv);
   if (!mKF->PropagateToHit(&trCorK1, &mKHit, kFALSE, kFALSE)) {
      return false;
   }
   trCorK1.SetDirection(MpdKalmanTrack::kInward);
   TVector3 m1 = trCorK1.Momentum3();

   hitTmp.SetPos(trCorK2.GetPos());
   trCorK2.SetParamNew(*trCorK2.GetParam());
   trCorK2.SetPos(trCorK2.GetPosNew());
   trCorK2.ReSetWeight();
   TMatrixDSym w = *trCorK1.GetWeight(); // save current weight matrix
   mKHit.SetPos(rConv);
   if (!mKF->PropagateToHit(&trCorK2, &mKHit, kFALSE, kFALSE)) {
      return false;
   }
   trCorK2.SetDirection(MpdKalmanTrack::kInward);
   TVector3 m2 = trCorK2.Momentum3();

   float qt, alpha;
   ArmenterosPodolanski(m1, m2, qt, alpha);

   mhArmPo->Fill(alpha, qt);
   if (isTrue) {
      mhArmPoTrue->Fill(alpha, qt);
   }
   if (!ArmenterosQtCut(qt, alpha, gamEE)) {
      // return false;
   }
   mhCutEff->Fill(17., pt);

   // Asymmetry cut
   float asym1 = m1.Mag() / gamEE.Momentum();
   float asym2 = m2.Mag() / gamEE.Momentum();
   mhAsym->Fill(asym1, pt);
   mhAsym->Fill(asym2, pt);
   if (isTrue) {
      mhAsymTrue->Fill(asym1, pt);
      mhAsymTrue->Fill(asym2, pt);
   }
   if (!(AsymmetryCut(asym1, pt) && AsymmetryCut(asym2, pt))) {
      return kFALSE;
   }
   mhCutEff->Fill(18., pt);

   float cospsi = CosPsiPair(m1, m2);
   mhCosPsi->Fill(cospsi, pt);
   if (isTrue) { // same for true electrontracks
      mhCosPsiTrue->Fill(cospsi, pt);
   }
   if (cospsi < mParams.mCosPsiCut) {
      return kFALSE;
   }
   mhCutEff->Fill(19., pt);

   // if(TMath::Abs(photonAOD->GetDCAzToPrimVtx()) > fDCAZPrimVtxCut) { //DCA Z cut of photon to primary vertex
   //   return kFALSE;
   // }
   mhCutEff->Fill(20., pt);

   // if(fHistoInvMassafter)fHistoInvMassafter->Fill(photon->GetMass());
   // if(fHistoArmenterosafter)fHistoArmenterosafter->Fill(photon->GetArmenterosAlpha(),photon->GetArmenterosQt());
   // if(fHistoPsiPairDeltaPhiafter)fHistoPsiPairDeltaPhiafter->Fill(deltaPhi,photon->GetPsiPair());
   // if(fHistoKappaafter)fHistoKappaafter->Fill(photon->GetPhotonPt(), GetKappaTPC(photon, event));
   // if(fHistoAsymmetryafter){
   //   if(photon->GetPhotonP()!=0 &&
   //   electronCandidate->P()!=0)fHistoAsymmetryafter->Fill(photon->GetPhotonP(),electronCandidate->P()/photon->GetPhotonP());
   // }

   mStorage.emplace_back(isTrue, pt, chi2, gamEE.GetMass(), rConv, angle, dist, qt, alpha, asym1, cospsi);

   v.SetXYZT((gamEE.Pt()) * TMath::Cos(gamEE.Phi()), (gamEE.Pt()) * TMath::Sin(gamEE.Phi()),
             TMath::Sign(TMath::Sqrt(gamEE.Momentum() * gamEE.Momentum() - gamEE.Pt() * gamEE.Pt()),
                         TMath::Cos(gamEE.Theta())),
             gamEE.Momentum());
   v.setPrimary(commonParentId);

   mhConvSp->Fill(v.Pt(), v.Eta());
   if (isTrue) {
      mhConvSpTrue->Fill(v.Pt(), v.Eta());
   }

   return true;
}
///________________________________________________________________________
bool MpdConvPi0::TestHybrid(MpdPhoton &c, MpdPhoton &v0) const
{
   double dphi = 999., dz = 999.;
   // Test if cluster match with any track from V0
   MpdEmcClusterKI *clu  = (MpdEmcClusterKI *)mEMCClusters->At(c.getTr1());
   double           xEMC = clu->GetX();
   double           yEMC = clu->GetY();
   double           zEMC = clu->GetZ();

   // int itr1 = v0.getTr1() ;
   // int itr2 = v0.getTr2() ;

   MpdTpcKalmanTrack *tr1 = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(v0.getTr1());
   MpdTpcKalmanTrack  tr1tmp(*tr1);
   tr1tmp.SetParam(*tr1tmp.GetParamAtHit());
   tr1tmp.SetParamNew(*tr1tmp.GetParamAtHit());
   tr1tmp.SetWeight(*tr1tmp.GetWeightAtHit());
   tr1tmp.SetPos(tr1tmp.GetPosAtHit());
   tr1tmp.SetPosNew(tr1tmp.GetPos());
   tr1tmp.SetLength(tr1tmp.GetLengAtHit());

   // Propagate to EMC cluser radius
   dphi = 999.;
   dz   = 999.;
   MpdKalmanHit hEnd;
   hEnd.SetType(MpdKalmanHit::kFixedR);
   double rClu = 170; // MpdEmcGeoUtils::GetInstance()->Rperp(zEMC) ;
   hEnd.SetPos(rClu);
   MpdKalmanFilter *pKF = MpdKalmanFilter::Instance("KF", "KF");

   if (pKF->PropagateToHit(&tr1tmp, &hEnd, kTRUE)) {

      double phi = tr1tmp.GetParamNew(0) / tr1tmp.GetPosNew();
      double z   = tr1tmp.GetParamNew(1);
      double r   = tr1tmp.GetPosNew();
      double x   = r * TMath::Cos(phi);
      double y   = r * TMath::Sin(phi);

      dphi = TMath::Sqrt((x - xEMC) * (x - xEMC) + (y - yEMC) * (y - yEMC));
      dz   = TMath::Abs(z - zEMC);
   }

   MpdTpcKalmanTrack *tr2 = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(v0.getTr2());
   MpdTpcKalmanTrack  tr2tmp(*tr2);
   tr2tmp.SetParam(*tr2tmp.GetParamAtHit());
   tr2tmp.SetParamNew(*tr2tmp.GetParamAtHit());
   tr2tmp.SetWeight(*tr2tmp.GetWeightAtHit());
   tr2tmp.SetPos(tr2tmp.GetPosAtHit());
   tr2tmp.SetPosNew(tr2tmp.GetPos());
   tr2tmp.SetLength(tr2tmp.GetLengAtHit());

   if (pKF->PropagateToHit(&tr2tmp, &hEnd, kTRUE)) {

      double phi = tr2tmp.GetParamNew(0) / tr2tmp.GetPosNew();
      double z   = tr2tmp.GetParamNew(1);
      double r   = tr2tmp.GetPosNew();
      double x   = r * TMath::Cos(phi);
      double y   = r * TMath::Sin(phi);

      double ddphi = TMath::Sqrt((x - xEMC) * (x - xEMC) + (y - yEMC) * (y - yEMC));
      double ddz   = TMath::Abs(z - zEMC);
      if (ddphi * ddphi + ddz * ddz < dphi * dphi + dz * dz) {
         dphi = ddphi;
         dz   = ddz;
      }
   }

   return (dphi * dphi / (25. * 25.) + dz * dz / (9. * 9.)) > 1.;
}
long int MpdConvPi0::IsSameParent(long int prim1, long int prim2) const
{
   // Looks through parents and finds if there was commont pi0 among ancestors

   if (!isMC) return -1; // can not say anything

   while (prim1 != -1) {
      long int pr2 = prim2;

      while (pr2 != -1) {
         if (prim1 == pr2) {
            return prim1;
         }
         pr2 = (static_cast<MpdMCTrack *>(mMCTracks->At(pr2)))->GetMotherId();
      }
      prim1 = (static_cast<MpdMCTrack *>(mMCTracks->At(prim1)))->GetMotherId();
   }
   return -1;
}

MpdHelix MpdConvPi0::MakeHelix(const MpdKalmanTrack &tr) const
{
   float r   = tr.GetPosNew();
   float phi = tr.GetParam(0) / r;
   float x   = r * TMath::Cos(phi);
   float y   = r * TMath::Sin(phi);
   float dip = tr.GetParam(3);
   float cur = 0.3 * 0.01 * 5 / 10; // 5 kG
   cur *= TMath::Abs(tr.GetParam(4));
   TVector3 o(x, y, tr.GetParam(1));
   Int_t    h = (Int_t)TMath::Sign(1.1, tr.GetParam(4));
   MpdHelix helix(cur, dip, tr.GetParam(2) - TMath::PiOver2() * h, o, h);
   return helix;
}

float MpdConvPi0::dEdx_sigma(float dEdx, float mom) const
{
   // To be moved to centralized class
   if (mom < 0.05) mom = 0.05;

   float mean[7]  = {9.793192e+002, -5.234570e-003, -3.178321e+000, -4.987832e-002,
                    3.617478e-002, -1.021387e-001, 9.169614e+002};
   float width[5] = {-1.589388e+003, 1.834372e+003, 4.125626e-003, 9.376418e-001, -1.466546e-001};

   float mean_exp, width_exp;

   mean_exp =
      mean[0] / mom / mom *
         (mean[1] * log(mom * mom) - mean[2] * mom * mom - mean[3] * mom - mean[4] - mean[5] * mom * mom * mom) +
      mean[6];
   width_exp = width[0] + width[1] * pow(mom, width[2]) + width[3] / pow(mom - width[4], 3);

   return (dEdx - mean_exp) / width_exp;
}

float MpdConvPi0::Beta_sigma(float beta, float mom) const
{
   // To be moved to centralized class
   if (mom < 0.05) mom = 0.05;

   float mean[7]  = {3.150000e+003, -4.833115e-008, 1.688117e+000, -7.840445e-007,
                    3.034956e-007, -3.852622e-007, 5.318573e+003};
   float width[5] = {5.854497e-003, 6.078866e-003, -1.174312e-001, 3.039271e-006, -8.370411e-002};

   float mean_exp, width_exp;

   mean_exp =
      mean[0] / mom / mom *
         (mean[1] * log(mom * mom) - mean[2] * mom * mom - mean[3] * mom - mean[4] - mean[5] * mom * mom * mom) +
      mean[6] - 0.001;
   width_exp = width[0] + width[1] * pow(mom, width[2]) + width[3] / pow(mom - width[4], 4);

   return (beta - mean_exp) / width_exp;
}

void MpdConvPi0::ArmenterosPodolanski(TVector3 &m1, TVector3 &m2, float &qt, float &alpha) const
{

   alpha = 0., qt = 0.;

   TVector3 s = m1 + m2;

   float pn  = m1.Mag();
   float pln = m1.Dot(s);
   float plp = m2.Dot(s);

   if (pn == 0.0) return;
   alpha    = (plp - pln) / (plp + pln);
   float sm = s.Mag();
   if (sm > 0) {
      qt = m1.Cross(s).Mag() / sm;
   }
}

//________________________________________________________________________
bool MpdConvPi0::ArmenterosQtCut(float qt, float alpha, MpdParticle &part) const
{ // Armenteros Qt Cut
   // if(mParams.mDo2DQt){
   //   if(mParams.mDoQtGammaSelection==1){
   //     if (
   //     !(TMath::Power(photon->GetArmenterosAlpha()/mParams.mMaxPhotonAsymmetry,2)+TMath::Power(photon->GetArmenterosQt()/mParams.mQtMax,2)
   //     < 1) ){
   //       return false;
   //     }
   //   } else if(mParams.mDoQtGammaSelection==2){
   //     float qtMaxPtDep = mParams.mQtPtMax*photon->GetPhotonPt();
   //     if (qtMaxPtDep > mParams.mQtMax)
   //       qtMaxPtDep      = mParams.mQtMax;
   //     if (
   //     !(TMath::Power(photon->GetArmenterosAlpha()/mParams.mMaxPhotonAsymmetry,2)+TMath::Power(photon->GetArmenterosQt()/qtMaxPtDep,2)
   //     < 1) ){
   //       return false;
   //     }
   //   }
   // } else {
   //   if(mParams.mDoQtGammaSelection==1){
   //     if(photon->GetArmenterosQt()>mParams.mQtMax){
   //       return false;
   //     }
   //   } else if(mParams.mDoQtGammaSelection==2){
   //     Float_t qtMaxPtDep = mParams.mQtPtMax*photon->GetPhotonPt();
   //     if (qtMaxPtDep > mParams.mQtMax)
   //       qtMaxPtDep      = mParams.mQtMax;
   //     if(photon->GetArmenterosQt()>qtMaxPtDep){
   //       return false;
   //     }
   //   }
   // }
   return true;
}

///________________________________________________________________________
bool MpdConvPi0::AsymmetryCut(float asym, float pt) const
{
   // Cut on Energy Asymmetry

   // for(Int_t ii=0;ii<2;ii++){

   //   AliVTrack *track=GetTrack(event,photon->GetTrackLabel(ii));

   //   if(fDoPhotonPDependentAsymCut){
   //     float trackNegAsy=0;
   //     if (photon->GetPhotonP()!=0.){
   //         trackNegAsy= track->P()/photon->GetPhotonP();
   //     }

   //     if( trackNegAsy > fFAsymmetryCut->Eval(photon->GetPhotonP()) || trackNegAsy
   //     < 1.-fFAsymmetryCut->Eval(photon->GetPhotonP()) ){
   //       return kFALSE;
   //     }

   //   } else {
   //     if( track->P() > fMinPPhotonAsymmetryCut ){
   //       float trackNegAsy=0;
   //       if (photon->GetPhotonP()!=0.){
   //         trackNegAsy= track->P()/photon->GetPhotonP();
   //       }

   //       if( trackNegAsy<fMinPhotonAsymmetry ||trackNegAsy>(1.- fMinPhotonAsymmetry)){
   //         return kFALSE;
   //       }
   //     }
   //   }

   // }
   return true;
}
///________________________________________________________________________
float MpdConvPi0::CosPsiPair(TVector3 &p1, TVector3 &p2) const
{

   // float p1[3] = {tr1->GetPx(),tr1->GetPy(),tr1->GetPz()};
   // float p2[3] = {tr2->GetPx(),tr2->GetPy(),tr2->GetPz()};
   // float u[3] = {p1[0]+p2[0],p1[1]+p2[1],p1[2]+p2[2]};
   TVector3 u = p1 + p2;

   // float normp1 = sqrt( (p1[0]*p1[0]) + (p1[1]*p1[1]) + (p1[2]*p1[2]) );
   // float normp2 = sqrt( (p2[0]*p2[0]) + (p2[1]*p2[1]) + (p2[2]*p2[2]) );
   // float normu  = sqrt( (u[0]*u[0]) + (u[1]*u[1]) + (u[2]*u[2]) );

   // for(int i=3; i--;){
   //   p1[i] /= normp1;
   //   p2[i] /= normp2;
   //   u[i] /= normu;
   // }

   TVector3 v = p1.Cross(p2);
   TVector3 w = u.Cross(v);
   TVector3 z(0, 0, 1.);
   TVector3 wc = u.Cross(z);
   return wc.Angle(w);
}

float MpdConvPi0::distCPV(float dphi, float dz, float E) const
{

   float sigmaPhi = 3.66601 - 4.63964e-01 / E + 2.08779e-01 / E / E;
   float sigmaZ   = 2.58409 - 1.87502e-01 / E + 2.40143e-01 / E / E;
   dphi           = dphi / sigmaPhi;
   dz             = dz / sigmaZ;
   return sqrt(dphi * dphi + dz * dz);
}

float MpdConvPi0::lambdaCut(float l1, float l2, float E) const
{

   float longM  = 4.28333;
   float shortM = 1.88168 - 5.06456e-01 * exp(-E / 3.83640e-01);
   float longS  = 1.05616 - 2.12212e-01 * exp(-E / 5.46530e-01);
   float shortS = 7.58640e-01 - 3.97720e-01 * exp(-E / 3.18150e-01);
   float c      = -1.0 + 5.42460e-01 * exp(-E / 3.22982e-01);

   return (l1 - longM) * (l1 - longM) / (longS * longS * 2.) + (l2 - shortM) * (l2 - shortM) / (shortS * shortS * 2.) +
          c * (l1 - longM) * (l2 - shortM) / (longS * shortS * 2.);
}

float MpdConvPi0::tofCut(float time, float E) const
{
   // return distance of time from expected photon arraival in sigma (with sign)
   float sigma = 1.86166 * TMath::Exp(-E / 0.0259728) + 0.347552; // resolution in ns
   return time / sigma;
}

float MpdConvPi0::Nonlinearity(float oldE) const
{

   float x = TMath::Min(oldE, 2.5f);
   return 2.9411765 * oldE / (0.97630219 + 7.194380e-002 * x - 4.491255e-002 * x * x + 8.362250e-003 * x * x * x);
}