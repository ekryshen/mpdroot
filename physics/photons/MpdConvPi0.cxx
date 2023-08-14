#include <iostream>
#include <fstream> // std::ifstream

#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"

#include "MpdKalmanFilter.h"
#include "MpdMCTrack.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdEmcClusterKI.h"
#include "MpdParticle.h"
#include "MpdVertex.h"
#include "MpdEvent.h"
#include "MpdEmcGeoUtils.h"
#include "MpdV0.h"
#include "MpdV0Maker.h"

#include "MpdConvPi0.h"

ClassImp(MpdConvPi0);

MpdConvPi0::MpdConvPi0(const char *name, const char *outputName) : MpdAnalysisTask(name, outputName)
{
   mParamConfig = Form("%s", outputName);
}

void MpdConvPi0::UserInit()
{
   cout << "MpdV0Maker::UserInit():\n";
   mParams.ReadFromFile(mParamConfig);
   mParams.Print();
   applySelection = mParams.mApplySelection;

   // Prepare histograms etc.
   fOutputList = new TList();
   fOutputList->SetOwner(kTRUE);

   mHistoCentBins = mCentBins.size() - 1;
   mAxCent        = new TAxis(mHistoCentBins, mCentBins.data());

   TH1::AddDirectory(false); // sets a global switch disabling the reference to histos in gROOT and their overwriting

   // General QA
   mhEventCutEff = addHist(
      new TH1F("hEventCutEff", "Number of events;Cut ID;N_{events}", eventCutId::kNcuts, 0, eventCutId::kNcuts));
   mhVertex     = addHist(new TH1F("hVertex", "Event vertex distribution;Z_{vtx} (cm);N_{events}", 50, -100., 100.));
   mhCentrality = addHist(new TH1F("hCentrality", "Centrality distribution;Centrality (%);N_{events}", 100, 0., 100.));

   // track selection
   const int    nPtbin = 100;
   const double pTmin  = 0.;
   const double pTmax  = 5.;
   mhTrackCutEff       = addHist(new TH2F("hTrackCutEff", "track cut efficiency;Cut ID;p_{T} (GeV/#it{c});Efficiency",
                                          trackCutId::kNcuts, 0., trackCutId::kNcuts, nPtbin, pTmin, pTmax));
   mhTrackNhits  = addHist(new TH1F("hTrackNhits", "Number of hits per track;N_{hits};N_{tracks}", 100, 0., 100.));
   mhTrackEtaPt  = addHist(new TH2F("hTrackEtaPt", "track occupancy pt vs eta;p_{T} (GeV/#it{c});#eta;N_{tracks}",
                                    nPtbin, pTmin, pTmax, 100, -1.5, 1.5));
   mhTrackProbEl = addHist(
      new TH2F("hTrackProbEl", "Electron probability;p (e^{-});p_{T} (GeV/#it{c})", 100, 0., 1., nPtbin, pTmin, pTmax));
   mhTrackNsigDEdx = addHist(new TH2F("hTrackNsigDEdx", "N#sigma(dE/dx);N#sigma(dE/dx);p_{T} (GeV/#it{c})", 100, -10,
                                      10., nPtbin, pTmin, pTmax));
   mhTrackNsigBeta = addHist(new TH2F("hTrackNsigBeta", "N#sigma(#beta);N#sigma(#beta);p_{T} (GeV/#it{c})", 100, -10,
                                      10., nPtbin, pTmin, pTmax));
   mhTrackNsigDEdxNsigBeta = addHist(
      new TH2F("hTrackNsigDEdxNsigBeta", "N#sigma(beta):N#sigma(dE/dx);N#sigma(dE/dx);N#sigma(beta);p_{T} (GeV/#it{c})",
               100, -10, 10., 100, -10, 10.));
   if (isMC) {
      mhTrackNhitsTrue            = addHistClone(mhTrackNhits, "True");
      mhTrackEtaPtTrue            = addHistClone(mhTrackEtaPt, "True");
      mhTrackProbElTrue           = addHistClone(mhTrackProbEl, "True");
      mhTrackNsigDEdxTrue         = addHistClone(mhTrackNsigDEdx, "True");
      mhTrackNsigBetaTrue         = addHistClone(mhTrackNsigBeta, "True");
      mhTrackNsigDEdxNsigBetaTrue = addHistClone(mhTrackNsigDEdxNsigBeta, "True");
   }

   // V0 selection
   mhV0CutEff = addHist(new TH2F("hV0CutEff", "V0 cut efficiency;Cut ID;p_{T} (GeV/#it{c});Efficiency", V0CutId::kNcuts,
                                 0., V0CutId::kNcuts, nPtbin, pTmin, pTmax));
   mhAlpha    = addHist(
         new TH2F("hAlpha", "#alpha distribution;#alpha (rad);p_{T} (GeV/#it{c})", 100, 0., 3.14, nPtbin, pTmin, pTmax));
   mhChi2    = addHist(new TH2F("hChi2", "#chi^{2};#chi^{2};p_{T} (GeV/#it{c})", 100, 0., 10., nPtbin, pTmin, pTmax));
   mhConvMap = addHist(new TH3F("hConvMap", "Conversion map (r,phi,z);r (cm);#phi (rad);z (cm)", 100, 0., 280., 100, 0.,
                                TMath::Pi(), 100, -200., 200.));
   mhV0rConv =
      addHist(new TH2F("hV0rConv", "R^{conv};R^{conv} (cm);p_{T} (GeV/#it{c})", 200, 0., 200., nPtbin, pTmin, pTmax));
   mhDist   = addHist(new TH2F("hDist", "track DCA;DCA (cm);p_{T} (GeV/#it{c})", 100, 0., 10., nPtbin, pTmin, pTmax));
   mhMassEE = addHist(
      new TH2F("hMassEE", "m_{ee};m_{ee} (GeV/#it{c}^{2});p_{T} (GeV/#it{c})", 100, 0., 0.3, nPtbin, pTmin, pTmax));
   mhCosPsi = addHist(new TH2F("cosPsi", "cos(#psi);cos(#psi);p_{T} (GeV/#it{c})", 100, -1., 1., nPtbin, pTmin, pTmax));
   mhArmPo  = addHist(new TH2F("Armenteros", "Armenteros", 100, -1, 1, 100, 0, 0.3));
   mhAsym   = addHist(new TH2F("Asymetry", "Asymmetry;Asymmetry;p_{T} (GeV/#it{c})", 200, 0, 1, nPtbin, pTmin, pTmax));
   mhConvSp = addHist(new TH2F("ConvSp", "Conv Sp;p_{T} (GeV/#it{c});#eta", nPtbin, pTmin, pTmax, 100, -1.5, 1.5));
   if (isMC) {
      mhConvMapTrue = addHistClone(mhConvMap, "True");
      mhAlphaTrue   = addHistClone(mhAlpha, "True");
      mhChi2True    = addHistClone(mhChi2, "True");
      mhV0rConvTrue = addHistClone(mhV0rConv, "True");
      mhDistTrue    = addHistClone(mhDist, "True");
      mhMassEETrue  = addHistClone(mhMassEE, "True");
      mhCosPsiTrue  = addHistClone(mhCosPsi, "True");
      mhArmPoTrue   = addHistClone(mhArmPo, "True");
      mhAsymTrue    = addHistClone(mhAsym, "True");
      mhConvSpTrue  = addHistClone(mhConvSp, "True");
   }

   // Cluster selection
   int   nBinsE = 100;
   float Emin = 0, Emax = 2.5;
   mhCluCutEff = addHist(new TH2F("hCluCutEff", "Cluster cut efficiency;Cut ID;e (GeV)", clustCutId::kNcuts, 0.,
                                  clustCutId::kNcuts, nBinsE, Emin, Emax));
   mhCluMult   = addHist(new TH2F("hCluMult", "Cluster multiplicity;Mult;e (GeV)", 100, 0, 100, nBinsE, Emin, Emax));
   mhCluTofSigma =
      addHist(new TH2F("hCluTofSigma", "Cluster #sigma_{TOF};#sigma_{TOF};e (GeV)", 200, -10, 10, nBinsE, Emin, Emax));
   mhCluDistCPV =
      addHist(new TH2F("hCluDistCPV", "Cluster DistCPV;DistCPV (cm);e (GeV)", 100, 0, 30, nBinsE, Emin, Emax));
   // Inv mass histos
   const int   nMbins = 200;
   const float mMax   = 1.;
   for (int cen = 0; cen < mHistoCentBins; cen++) {
      std::string centrality = Form("%.0f-%.0f", mAxCent->GetBinLowEdge(cen + 1), mAxCent->GetBinUpEdge(cen + 1));
      mhRealCalo.push_back(addHist(
         new TH2F(Form("hRealCalo_cen_%s", centrality.c_str()),
                  Form("Real inv mass, calorimeter, %s;M_{inv} (GeV/#it{c}));p_{T} (GeV/#it{c})", centrality.c_str()),
                  nMbins, 0., mMax, nPtbin, pTmin, pTmax)));
      mhMixedCalo.push_back(addHist(
         new TH2F(Form("hMixedCalo_cen_%s", centrality.c_str()),
                  Form("Mixed inv mass, calorimeter, %s;M_{inv} (GeV/#it{c});p_{T} (GeV/#it{c})", centrality.c_str()),
                  nMbins, 0., mMax, nPtbin, pTmin, pTmax)));
      mhRealHybrid.push_back(
         addHist(new TH2F(Form("hRealHybrid_cen_%s", centrality.c_str()),
                          Form("Real inv mass, hybrid, %s;M_{inv} (GeV/#it{c});p_{T} (GeV/#it{c})", centrality.c_str()),
                          nMbins, 0., mMax, nPtbin, pTmin, pTmax)));
      mhMixedHybrid.push_back(addHist(
         new TH2F(Form("hMixedHybrid_cen_%s", centrality.c_str()),
                  Form("Mixed inv mass, hybrid, %s;M_{inv} (GeV/#it{c});p_{T} (GeV/#it{c})", centrality.c_str()),
                  nMbins, 0., mMax, nPtbin, pTmin, pTmax)));
      mhRealConv.push_back(addHist(
         new TH2F(Form("hRealConv_cen_%s", centrality.c_str()),
                  Form("Real inv mass, conversion, %s;M_{inv} (GeV/#it{c});p_{T} (GeV/#it{c})", centrality.c_str()),
                  nMbins, 0., mMax, nPtbin, pTmin, pTmax)));
      mhMixedConv.push_back(addHist(
         new TH2F(Form("hMixedConv_cen_%s", centrality.c_str()),
                  Form("Mixed inv mass, conversion, %s;M_{inv} (GeV/#it{c});p_{T} (GeV/#it{c})", centrality.c_str()),
                  nMbins, 0., mMax, nPtbin, pTmin, pTmax)));

      // flow
      mp3V1etaPtMinvCalo.push_back(
         addHist(new TProfile3D(Form("p3V1etaPtMinvCalo_cen_%s", centrality.c_str()),
                                Form("V1, Calo, %s;#eta;p_{T} (GeV/#it{c});M_{inv} (GeV/#it{c}))", centrality.c_str()),
                                15, -1.5, 1.5, 25, 0, 5, nMbins, 0., mMax)));
      mp3V1etaPtMinvHybrid.push_back(addHist(
         new TProfile3D(Form("p3V1etaPtMinvHybrid_cen_%s", centrality.c_str()),
                        Form("V1, Hybrid, %s;#eta;p_{T} (GeV/#it{c});M_{inv} (GeV/#it{c}))", centrality.c_str()), 15,
                        -1.5, 1.5, 25, 0, 5, nMbins, 0., mMax)));
      mp3V1etaPtMinvConv.push_back(
         addHist(new TProfile3D(Form("p3V1etaPtMinvConv_cen_%s", centrality.c_str()),
                                Form("V1, Conv, %s;#eta;p_{T} (GeV/#it{c});M_{inv} (GeV/#it{c}))", centrality.c_str()),
                                15, -1.5, 1.5, 25, 0, 5, nMbins, 0., mMax)));

      if (isMC) {
         mhRealCaloTruePi.push_back(addHistClone(mhRealCalo.back(), "TruePi"));
         mhRealHybridTruePi.push_back(addHistClone(mhRealHybrid.back(), "TruePi"));
         mhRealConvTruePi.push_back(addHistClone(mhRealConv.back(), "TruePi"));
         mhRealCaloTruePh.push_back(addHistClone(mhRealCalo.back(), "TruePh"));
         mhRealHybridTruePh.push_back(addHistClone(mhRealHybrid.back(), "TruePh"));
         mhRealConvTruePh.push_back(addHistClone(mhRealConv.back(), "TruePh"));
         mhRealCaloTrueAll.push_back(addHistClone(mhRealCalo.back(), "TrueAll"));
         mhRealHybridTrueAll.push_back(addHistClone(mhRealHybrid.back(), "TrueAll"));
         mhRealConvTrueAll.push_back(addHistClone(mhRealConv.back(), "TrueAll"));
         hPrimPi.push_back(
            addHist(new TH2F(Form("PrimaryPi0_cen_%s", centrality.c_str()),
                             Form("Primary Pi0, centrality %s;#eta;p_{T} (GeV/#it{c})", centrality.c_str()), 100, -2.,
                             2., nPtbin, pTmin, pTmax)));
         hPrimPh.push_back(
            addHist(new TH2F(Form("PrimaryPhoton_cen_%s", centrality.c_str()),
                             Form("Primary Photon, centrality %s;#eta;p_{T} (GeV/#it{c})", centrality.c_str()), 100,
                             -2., 2., nPtbin, pTmin, pTmax)));
         mp3V1etaPtMinvCaloTruePi.push_back(addHistClone(mp3V1etaPtMinvCalo.back(), "TruePi"));
         mp3V1etaPtMinvHybridTruePi.push_back(addHistClone(mp3V1etaPtMinvHybrid.back(), "TruePi"));
         mp3V1etaPtMinvConvTruePi.push_back(addHistClone(mp3V1etaPtMinvConv.back(), "TruePi"));
         mp3V1etaPtMinvCaloTruePh.push_back(addHistClone(mp3V1etaPtMinvCalo.back(), "TruePh"));
         mp3V1etaPtMinvHybridTruePh.push_back(addHistClone(mp3V1etaPtMinvHybrid.back(), "TruePh"));
         mp3V1etaPtMinvConvTruePh.push_back(addHistClone(mp3V1etaPtMinvConv.back(), "TruePh"));
         mp2V1etaPtPrimPiEP.push_back(addHist(
            new TProfile2D(Form("p2V1etaPt_cen_%s_PiEP", centrality.c_str()),
                           Form("V1 Pi EP, True, %s;#eta;p_{T} (GeV/#it{c});M_{inv} (GeV/#it{c}))", centrality.c_str()),
                           15, -1.5, 1.5, 25, 0, 5)));
         mp2V1etaPtPrimPhEP.push_back(addHist(
            new TProfile2D(Form("p2V1etaPt_cen_%s_PhEP", centrality.c_str()),
                           Form("V1 Ph EP, True, %s;#eta;p_{T} (GeV/#it{c});M_{inv} (GeV/#it{c}))", centrality.c_str()),
                           15, -1.5, 1.5, 25, 0, 5)));
         mp2V1etaPtPrimPiRP.push_back(addHist(
            new TProfile2D(Form("p2V1etaPt_cen_%s_PiRP", centrality.c_str()),
                           Form("V1 Pi RP, True, %s;#eta;p_{T} (GeV/#it{c});M_{inv} (GeV/#it{c}))", centrality.c_str()),
                           15, -1.5, 1.5, 25, 0, 5)));
         mp2V1etaPtPrimPhRP.push_back(addHist(
            new TProfile2D(Form("p2V1etaPt_cen_%s_PhRP", centrality.c_str()),
                           Form("V1 Ph RP, True, %s;#eta;p_{T} (GeV/#it{c});M_{inv} (GeV/#it{c}))", centrality.c_str()),
                           15, -1.5, 1.5, 25, 0, 5)));
      }
   }

   cout << "MpdV0Maker::UserInit(): complete\n";
}
//--------------------------------------
void MpdConvPi0::ProcessEvent(MpdAnalysisEvent &event)
{
   mMpdGlobalTracks = event.fMPDEvent->GetGlobalTracks();
   mKalmanTracks    = event.fTPCKalmanTrack;
   mEMCClusters     = event.fEMCCluster;
   mMCTracks        = event.fMCTrack;
   psiRP            = event.fMCEventHeader->GetRotZ();

   if (applySelection && !selectEvent(event)) {
      return;
   }

   int nTracks = mMpdGlobalTracks->GetEntriesFast();
   isGoodTrack.resize(nTracks, true);
   for (int i = 0; i < nTracks; i++) {
      auto *track = static_cast<MpdTrack *>(mMpdGlobalTracks->At(i));
      if (!selectTrack(track)) isGoodTrack.at(i) = false;
   }

   if (isMC) {
      for (int i = 0; i < mMCTracks->GetEntriesFast(); i++) {
         auto *p = static_cast<MpdMCTrack *>(mMCTracks->At(i));
         if (p->GetStartX() * p->GetStartX() + p->GetStartY() * p->GetStartY() <
             1.) { // make configurable??? what about StartZ?
            TVector3 mom;
            p->GetMomentum(mom);
            float eta = mom.Eta(), pt = mom.Pt(), v1ep = cos(mom.Phi() - psiEP), v1rp = cos(mom.Phi() - psiRP);
            if (p->GetPdgCode() == 111) {
               hPrimPi.at(mCenBin)->Fill(eta, pt);
               mp2V1etaPtPrimPiEP.at(mCenBin)->Fill(eta, pt, v1ep);
               mp2V1etaPtPrimPiRP.at(mCenBin)->Fill(eta, pt, v1rp);
            }
            if (p->GetPdgCode() == 22) {
               hPrimPh.at(mCenBin)->Fill(eta, pt);
               mp2V1etaPtPrimPhEP.at(mCenBin)->Fill(eta, pt, v1ep);
               mp2V1etaPtPrimPhRP.at(mCenBin)->Fill(eta, pt, v1rp);
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

   mhEventCutEff->Fill(eventCutId::kNone);
   // first test if event filled?
   if (!event.fVertex) return false;
   mhEventCutEff->Fill(eventCutId::kVertexPresent);

   // Vertex z coordinate
   MpdVertex *vertex = (MpdVertex *)event.fVertex->First();
   vertex->Position(mPrimaryVertex);
   mhVertex->Fill(mPrimaryVertex.Z());
   if (applySelection && fabs(mPrimaryVertex.Z()) > mParams.mZvtxCut) return false;
   mhEventCutEff->Fill(eventCutId::kVertexZ);

   // Centrality
   float centrality = event.getCentrTPC();
   mhCentrality->Fill(centrality);
   mCenBin = mAxCent->FindBin(centrality) - 1;
   if (applySelection && !(mCenBin > 0 && mCenBin < mAxCent->GetNbins())) return false;
   mhEventCutEff->Fill(eventCutId::kGoodCentrality);

   // vtxZ
   mZvtxBin = 0.5 * nMixEventZ * (1 + mPrimaryVertex.Z() / mParams.mZvtxCut);
   if (mZvtxBin < 0) mZvtxBin = 0;
   if (mZvtxBin >= nMixEventZ) mZvtxBin = nMixEventZ - 1;

   // EP
   psiEP  = event.fMpdEP.GetPhiEP_FHCal_F_all();
   psiEPn = event.fMpdEP.GetPhiEP_FHCal_N_all();
   psiEPs = event.fMpdEP.GetPhiEP_FHCal_S_all();
   mEPBin = 0.5 * nMixEventEP * (1 + psiEP / TMath::Pi());

   return true;
}
//--------------------------------------
void MpdConvPi0::selectConversion(MpdAnalysisEvent &event)
{
   // Select V0 in this event
   mV0.clear();

   auto v0array = event.fV0;
   if (!v0array) {
      cout << "WARNING: MpdConvPi0::selectConversion: empty V0 array!\n";
      return;
   }
   long int nv0 = v0array->GetEntriesFast();
   for (long int i = 0; i < nv0; i++) {
      auto v0 = (MpdV0 *)v0array->UncheckedAt(i);
      if (selectV0(v0)) {
         mV0.emplace_back(v0->Px(), v0->Py(), v0->Pz(), v0->E());
         mV0.back().setPrimary(v0->getCommonParent());
         mV0.back().setTr1(v0->getTr1());
         mV0.back().setTr2(v0->getTr2());
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
      mhCluCutEff->Fill(clustCutId::kNone, e);

      if (e < mParams.mCluEmin) {
         continue;
      }
      mhCluCutEff->Fill(clustCutId::kE, e);

      int cluMult = clu->GetMultiplicity();
      mhCluMult->Fill(cluMult, e);
      if (cluMult < mParams.mCluMult) {
         continue;
      }
      mhCluCutEff->Fill(clustCutId::kMult, e);

      double dx   = clu->GetX() - mPrimaryVertex.X();
      double dy   = clu->GetY() - mPrimaryVertex.Y();
      double dz   = clu->GetZ() - mPrimaryVertex.Z();
      double r    = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
      double time = clu->GetTime() - r / 29979245800. * 1.e+9; // Time in ns

      float cluTofSigma = tofCut(time, e);
      mhCluTofSigma->Fill(cluTofSigma, e);
      if (fabs(cluTofSigma) > mParams.mCluTof) {
         continue;
      }
      mhCluCutEff->Fill(clustCutId::kTof, e);

      float l1, l2; // Dispersion axis
      clu->GetLambdas(l1, l2);
      // correct for z
      float izMax = TMath::Abs(32. * clu->GetZ());
      l1          = l1 * 1.6 / (1.6 + 0.0002 * izMax * izMax);
      l2          = l2 * 3.2 / (3.2 + 0.000023 * izMax * izMax * izMax);

      //    if(e>mParams.mCluDispEmin &&  lambdaCut(l1,l2,e) > mParams.mCluDisp){
      //      continue ;
      //    }
      mhCluCutEff->Fill(clustCutId::kDisp, e);

      float pp0      = par0[0] + par0[1] * mPrimaryVertex.Z();
      float pp1      = par1[0] + par1[1] * mPrimaryVertex.Z();
      float z_shift1 = pp0 + pp1 * log(e);

      float cluDistCPV = distCPV(clu->GetDPhi(), clu->GetDZ() + z_shift1, e);
      mhCluDistCPV->Fill(cluDistCPV, e);
      if (cluDistCPV < mParams.mCluCPV) {
         continue;
      }
      mhCluCutEff->Fill(clustCutId::kCPV, e);

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
         float          eta = sum.Eta(), pt = sum.Pt(), m = sum.M(), v1 = (sum.Phi() - psiEP);
         mhRealCalo.at(mCenBin)->Fill(m, pt);
         mp3V1etaPtMinvCalo.at(mCenBin)->Fill(eta, pt, m, v1);
         long int ip = MpdV0Maker::FindCommonParent(mClusters[i].primary(), mClusters[j].primary(), mMCTracks);
         if (ip >= 0) {
            int pdg = abs(static_cast<MpdMCTrack *>(mMCTracks->At(ip))->GetPdgCode());
            mhRealCaloTrueAll.at(mCenBin)->Fill(m, pt);
            if (pdg == 111) {
               mp3V1etaPtMinvCaloTruePi.at(mCenBin)->Fill(eta, pt, m, v1);
               mhRealCaloTruePi.at(mCenBin)->Fill(m, pt);
            }
            if (pdg == 22) {
               mp3V1etaPtMinvCaloTruePh.at(mCenBin)->Fill(eta, pt, m, v1);
               mhRealCaloTruePh.at(mCenBin)->Fill(m, pt);
            }
         }
      }
   }

   for (int i = 0; i < nClu; i++) {
      for (int j = 0; j < nV0; j++) {
         if (!TestHybrid(mClusters[i], mV0[j])) continue;
         TLorentzVector sum = mClusters[i] + mV0[j];
         float          eta = sum.Eta(), pt = sum.Pt(), m = sum.M(), v1 = (sum.Phi() - psiEP);
         mhRealHybrid.at(mCenBin)->Fill(m, pt);
         mp3V1etaPtMinvHybrid.at(mCenBin)->Fill(eta, pt, m, v1);
         long int ip = MpdV0Maker::FindCommonParent(mClusters[i].primary(), mV0[j].primary(), mMCTracks);
         if (ip >= 0) {
            int pdg = abs(static_cast<MpdMCTrack *>(mMCTracks->At(ip))->GetPdgCode());
            mhRealHybridTrueAll.at(mCenBin)->Fill(m, pt);
            if (pdg == 111) {
               mp3V1etaPtMinvHybridTruePi.at(mCenBin)->Fill(eta, pt, m, v1);
               mhRealHybridTruePi.at(mCenBin)->Fill(m, pt);
            }
            if (pdg == 22) {
               mp3V1etaPtMinvHybridTruePh.at(mCenBin)->Fill(eta, pt, m, v1);
               mhRealHybridTruePh.at(mCenBin)->Fill(m, pt);
            }
         }
      }
   }

   for (int i = 0; i < nV0 - 1; i++) {
      for (int j = i + 1; j < nV0; j++) {
         TLorentzVector sum = mV0[i] + mV0[j];
         float          eta = sum.Eta(), pt = sum.Pt(), m = sum.M(), v1 = (sum.Phi() - psiEP);
         mhRealConv.at(mCenBin)->Fill(m, pt);
         mp3V1etaPtMinvConv.at(mCenBin)->Fill(eta, pt, m, v1);
         long int ip = MpdV0Maker::FindCommonParent(mV0[i].primary(), mV0[j].primary(), mMCTracks);
         if (ip >= 0) {
            int pdg = abs(static_cast<MpdMCTrack *>(mMCTracks->At(ip))->GetPdgCode());
            mhRealConvTrueAll.at(mCenBin)->Fill(m, pt);
            if (pdg == 111) {
               mp3V1etaPtMinvConvTruePi.at(mCenBin)->Fill(eta, pt, m, v1);
               mhRealConvTruePi.at(mCenBin)->Fill(m, pt);
            }
            if (pdg == 22) {
               mp3V1etaPtMinvConvTruePh.at(mCenBin)->Fill(eta, pt, m, v1);
               mhRealConvTruePh.at(mCenBin)->Fill(m, pt);
            }
         }
      }
   }

   // Mixed
   // calculate bin from zVertex-centrality-reaction plane
   int mixBin = mZvtxBin * nMixEventCent * nMixEventEP + mCenBin * nMixEventEP + mEPBin;

   for (auto &vm : mMixClu.at(mixBin)) {
      for (auto &v : mClusters) {
         TLorentzVector sum = v + vm;
         mhMixedCalo.at(mCenBin)->Fill(sum.M(), sum.Pt());
      }
   }
   for (auto &vm : mMixV0.at(mixBin)) {
      for (auto &v : mClusters) {
         TLorentzVector sum = v + vm;
         mhMixedHybrid.at(mCenBin)->Fill(sum.M(), sum.Pt());
      }
   }
   for (auto &vm : mMixClu.at(mixBin)) {
      for (auto &v : mV0) {
         TLorentzVector sum = v + vm;
         mhMixedHybrid.at(mCenBin)->Fill(sum.M(), sum.Pt());
      }
   }
   for (auto &vm : mMixV0.at(mixBin)) {
      for (auto &v : mV0) {
         TLorentzVector sum = v + vm;
         mhMixedConv.at(mCenBin)->Fill(sum.M(), sum.Pt());
      }
   }

   // Append new particles to queue and remove those at the beginning
   for (auto &v : mV0) {
      mMixV0.at(mixBin).emplace_back(v);
   }
   while (mMixV0.at(mixBin).size() > (UInt_t)mMaxMixSize) {
      mMixV0.at(mixBin).pop_front();
   }
   for (auto &v : mClusters) {
      mMixClu.at(mixBin).emplace_back(v);
   }
   while (mMixClu.at(mixBin).size() > (UInt_t)mMaxMixSize) {
      mMixClu.at(mixBin).pop_front();
   }
}

bool MpdConvPi0::selectTrack(MpdTrack *mpdtrack)
{
   bool isElectron = false;
   if (isMC) { // same for true electron tracks
      long int matchId = mpdtrack->GetID();
      if (matchId >= 0) {
         isElectron = abs((static_cast<MpdMCTrack *>(mMCTracks->At(matchId)))->GetPdgCode()) == 11;
      }
   }
   float pt = mpdtrack->GetPt(), eta = mpdtrack->GetEta();
   mhTrackCutEff->Fill(trackCutId::kNone, pt);

   int nHits = mpdtrack->GetNofHits();
   if (isElectron) mhTrackNhitsTrue->Fill(nHits, pt);
   if (applySelection && nHits < mParams.mNofHitsCut) return false;
   mhTrackNhits->Fill(nHits);
   mhTrackCutEff->Fill(trackCutId::kNhits, pt);

   if (isElectron) mhTrackEtaPtTrue->Fill(eta, pt);
   if (applySelection && !(fabs(eta) < mParams.mEtaCut && pt > mParams.mPtminCut)) return false;
   mhTrackEtaPt->Fill(eta, pt);
   mhTrackCutEff->Fill(trackCutId::kEtaPt, pt);

   float probEl = mpdtrack->GetPidProbElectron(); // TPC only or combined??? (actually none works)
   if (isElectron) mhTrackProbElTrue->Fill(probEl, pt);
   if (applySelection && probEl < mParams.mProbElCut) return false;
   mhTrackProbEl->Fill(probEl, fabs(mpdtrack->GetPt()));
   mhTrackCutEff->Fill(trackCutId::kProbElectron, pt);

   bool  hasTofHit   = (fabs(mpdtrack->GetTofDphiSigma()) < 3.0 && fabs(mpdtrack->GetTofDzSigma()) < 3.0);
   bool  tofElectron = false, dEdxElectron = false;
   float nSigma_dEdx = mpdtrack->GetTPCNSigma(kEl);
   float nSigma_beta = mpdtrack->GetTOFNSigma(kEl);
   if (isElectron) { // same for true electron tracks
      mhTrackNsigDEdxTrue->Fill(nSigma_dEdx, pt);
      mhTrackNsigBetaTrue->Fill(nSigma_beta, pt);
      mhTrackNsigDEdxNsigBetaTrue->Fill(nSigma_dEdx, nSigma_beta);
   }
   if (fabs(nSigma_dEdx) < mParams.mdEdxSigmaCut) {
      dEdxElectron = true;
   }
   if (hasTofHit && fabs(nSigma_beta) < mParams.mBetaSigmaCut) {
      tofElectron = true;
   }
   if (applySelection) {
      if (!dEdxElectron) return false;
      if (mParams.mRequireTOFpid && !tofElectron) return false;
   }
   mhTrackNsigDEdx->Fill(nSigma_dEdx, pt);
   mhTrackNsigBeta->Fill(nSigma_beta, pt);
   mhTrackNsigDEdxNsigBeta->Fill(nSigma_dEdx, nSigma_beta);
   mhTrackCutEff->Fill(trackCutId::kTpcTofPid, pt);

   return true;
}

bool MpdConvPi0::selectV0(MpdV0 *v0)
{
   // Check V0
   float pt = v0->Pt();
   mhV0CutEff->Fill(V0CutId::kNone, pt);

   if (applySelection && !(isGoodTrack.at(v0->getTr1()) && isGoodTrack.at(v0->getTr2()))) return false;
   mhV0CutEff->Fill(V0CutId::kGoodTracks, pt);

   if (applySelection && pt < 0.005) return false; // to avoid fpe
   mhV0CutEff->Fill(V0CutId::kNonZeroPt, pt);

   bool     isTrue         = false; // is true conv pair?
   long int commonParentId = -1;
   if (isMC) { // same for true electrontracks
      commonParentId = v0->getCommonParent();
      if (commonParentId >= 0) { // there is common parent
         isTrue = (static_cast<MpdMCTrack *>(mMCTracks->At(commonParentId))->GetPdgCode() == 22);
      }
   }

   float chi2 = v0->getChi2();
   if (isTrue) {
      mhChi2True->Fill(chi2, pt);
   }
   if (applySelection && chi2 > mParams.mChi2Cut) return false;
   mhChi2->Fill(chi2, pt);
   mhV0CutEff->Fill(V0CutId::kChi2, pt);

   float x, y, z;
   v0->getVertex(x, y, z);

   mhConvMap->Fill(x, y, z);
   if (isTrue) {
      mhConvMapTrue->Fill(x, y, z);
   }

   float rConv = v0->getRconv();
   if (isTrue) { // same for true electrontracks
      mhV0rConvTrue->Fill(rConv, pt);
   }
   if (applySelection && (rConv < mParams.mMinR2Cut || rConv > mParams.mMaxR2Cut)) return false;
   mhV0rConv->Fill(rConv, pt);
   mhV0CutEff->Fill(V0CutId::kRconv, pt);

   float alpha, qt;
   v0->getArmenteros(alpha, qt);
   if (isTrue) { // same for true electrontracks
      mhAlphaTrue->Fill(alpha, pt);
   }
   if (applySelection && alpha > mParams.mAlphaCut) return false;
   mhAlpha->Fill(alpha, pt);
   mhV0CutEff->Fill(V0CutId::kAlpha, pt);

   // if(applySelection &&  ePos->R() <= ((TMath::Abs(ePos->Vz()) * fLineCutZRSlope) - fLineCutZValue)){
   //   return false;  // line cut to exclude regions where we do not reconstruct
   // } else if (applySelection &&  fEtaCutMin != -0.1 &&   ePos->R() >= ((TMath::Abs(ePos->Vz()) * fLineCutZRSlopeMin)
   // - fLineCutZValueMin)){
   //   return false;
   // }

   float dist = v0->getDaughterDCA();

   if (isTrue) { // same for true electrontracks
      mhDistTrue->Fill(dist, pt);
   }
   if (applySelection && dist > mParams.mDistCut) return false;
   mhDist->Fill(dist, pt);
   mhV0CutEff->Fill(V0CutId::kDist, pt);

   float mass = v0->getMass();
   if (isTrue) { // same for true electrontracks
      mhMassEETrue->Fill(mass, pt);
   }
   if (applySelection && mass > mParams.mMassCut) return false;
   mhMassEE->Fill(mass, pt);
   mhV0CutEff->Fill(V0CutId::kMass, pt);

   if (isTrue) {
      mhArmPoTrue->Fill(alpha, qt);
   }
   if (applySelection && !ArmenterosQtCut(v0)) {
      // return false;
   }
   mhArmPo->Fill(alpha, qt);
   mhV0CutEff->Fill(V0CutId::kQt, pt);

   // Asymmetry cut
   float asym1, asym2;
   v0->getAsymmetry(asym1, asym2);
   if (isTrue) {
      mhAsymTrue->Fill(asym1, pt);
      mhAsymTrue->Fill(asym2, pt);
   }
   if (applySelection && !(AsymmetryCut(asym1, pt) && AsymmetryCut(asym2, pt))) {
      return false;
   }
   mhAsym->Fill(asym1, pt);
   mhAsym->Fill(asym2, pt);
   mhV0CutEff->Fill(V0CutId::kAsymmetry, pt);

   float cospsi = v0->getCospsi();
   if (isTrue) { // same for true electrontracks
      mhCosPsiTrue->Fill(cospsi, pt);
   }
   if (applySelection && cospsi < mParams.mCosPsiCut) return false;
   mhCosPsi->Fill(cospsi, pt);
   mhV0CutEff->Fill(V0CutId::kCosPsi, pt);

   // if(applySelection && TMath::Abs(photonAOD->GetDCAzToPrimVtx()) > fDCAZPrimVtxCut) { //DCA Z cut of photon to
   // primary vertex
   //   return false;
   // }
   mhV0CutEff->Fill(V0CutId::kPhotonAOD, pt);

   // if(fHistoInvMassafter)fHistoInvMassafter->Fill(photon->GetMass());
   // if(fHistoArmenterosafter)fHistoArmenterosafter->Fill(photon->GetArmenterosAlpha(),photon->GetArmenterosQt());
   // if(fHistoPsiPairDeltaPhiafter)fHistoPsiPairDeltaPhiafter->Fill(deltaPhi,photon->GetPsiPair());
   // if(fHistoKappaafter)fHistoKappaafter->Fill(photon->GetPhotonPt(), GetKappaTPC(photon, event));
   // if(fHistoAsymmetryafter){
   //   if(photon->GetPhotonP()!=0 &&
   //   electronCandidate->P()!=0)fHistoAsymmetryafter->Fill(photon->GetPhotonP(),electronCandidate->P()/photon->GetPhotonP());
   // }

   mhConvSp->Fill(v0->Pt(), v0->Eta());
   if (isTrue) {
      mhConvSpTrue->Fill(v0->Pt(), v0->Eta());
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

//________________________________________________________________________
bool MpdConvPi0::ArmenterosQtCut(MpdV0 *v0) const
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
   //       return false;
   //     }

   //   } else {
   //     if( track->P() > fMinPPhotonAsymmetryCut ){
   //       float trackNegAsy=0;
   //       if (photon->GetPhotonP()!=0.){
   //         trackNegAsy= track->P()/photon->GetPhotonP();
   //       }

   //       if( trackNegAsy<fMinPhotonAsymmetry ||trackNegAsy>(1.- fMinPhotonAsymmetry)){
   //         return false;
   //       }
   //     }
   //   }

   // }
   return true;
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
   // return distance of time from expected photon arrival in sigma (with sign)
   float sigma = 1.86166 * TMath::Exp(-E / 0.0259728) + 0.347552; // resolution in ns
   return time / sigma;
}

float MpdConvPi0::Nonlinearity(float oldE) const
{

   float x = TMath::Min(oldE, 2.5f);
   return 2.9411765 * oldE / (0.97630219 + 7.194380e-002 * x - 4.491255e-002 * x * x + 8.362250e-003 * x * x * x);
}
