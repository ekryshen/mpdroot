#include <iostream>
#include <fstream> // std::ifstream

#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TRandom.h"

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
   cout << "MpdConvPi0::UserInit():\n";
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
   hEventCutEff = addHist(
      new TH1F("hEventCutEff", "Number of events;Cut ID;N_{events}", eventCutId::kNcuts, 0, eventCutId::kNcuts));
   hVertex     = addHist(new TH1F("hVertex", "Event vertex distribution;Z_{vtx} (cm);N_{events}", 50, -100., 100.));
   hCentrality = addHist(new TH1F("hCentrality", "Centrality distribution;Centrality (%);N_{events}", 100, 0., 100.));

   // track selection
   const int    ptNbins = 100, etaNbins = 100;
   const double ptMin = 0., ptMax = 5., etaMin = -1.5, etaMax = 1.5;
   h3trackEtaPtCutEff =
      addHist(new TH3F("h3trackEtaPtCutEff", "track cut efficiency;#eta;p_{T} (GeV/#it{c});Cut ID", etaNbins, etaMin,
                       etaMax, ptNbins, ptMin, ptMax, trackCutId::kNcuts, 0., trackCutId::kNcuts));
   h3trackEtaPtNhits =
      addHist(new TH3F("h3trackEtaPtNhits", "Number of hits per track;#eta;p_{T} (GeV/#it{c});N_{hits};N_{tracks}",
                       etaNbins, etaMin, etaMax, ptNbins, ptMin, ptMax, 100, 0., 100.));
   if (isMC) h3trackEtaPtNhitsTrue = addHistClone(h3trackEtaPtNhits, "True");
   h2trackEtaPt = addHist(new TH2F("h2trackEtaPt", "track occupancy pt vs eta;#eta;p_{T} (GeV/#it{c});N_{tracks}",
                                   etaNbins, etaMin, etaMax, ptNbins, ptMin, ptMax));
   if (isMC) h2trackEtaPtTrue = addHistClone(h2trackEtaPt, "True");
   h3trackEtaPtProbEl = addHist(new TH3F("h3trackEtaPtProbEl", "Electron probability;#eta;p_{T} (GeV/#it{c});p (e^{-})",
                                         etaNbins, etaMin, etaMax, ptNbins, ptMin, ptMax, 100, 0., 1.));
   if (isMC) h3trackEtaPtProbElTrue = addHistClone(h3trackEtaPtProbEl, "True");
   h3trackEtaPtNsigDEdx =
      addHist(new TH3F("h3trackEtaPtNsigDEdx", "N#sigma(dE/dx);#eta;p_{T} (GeV/#it{c});N#sigma(dE/dx)", etaNbins,
                       etaMin, etaMax, ptNbins, ptMin, ptMax, 100, -10, 10.));
   if (isMC) h3trackEtaPtNsigDEdxTrue = addHistClone(h3trackEtaPtNsigDEdx, "True");
   h3trackEtaPtNsigBeta =
      addHist(new TH3F("h3trackEtaPtNsigBeta", "N#sigma(#beta);#eta;p_{T} (GeV/#it{c});N#sigma(#beta)", etaNbins,
                       etaMin, etaMax, ptNbins, ptMin, ptMax, 100, -10, 10.));
   if (isMC) h3trackEtaPtNsigBetaTrue = addHistClone(h3trackEtaPtNsigBeta, "True");
   h2trackNsigDEdxNsigBeta =
      addHist(new TH2F("h2trackNsigDEdxNsigBeta", "N#sigma(beta):N#sigma(dE/dx);N#sigma(dE/dx);N#sigma(beta)", 100, -10,
                       10., 100, -10, 10.));
   if (isMC) h2trackNsigDEdxNsigBetaTrue = addHistClone(h2trackNsigDEdxNsigBeta, "True");

   // V0 selection
   h3v0EtaPtCutEff = addHist(new TH3F("h3v0EtaPtCutEff", "V0 cut efficiency;#eta;p_{T} (GeV/#it{c});Cut ID", etaNbins,
                                      etaMin, etaMax, ptNbins, ptMin, ptMax, V0CutId::kNcuts, 0., V0CutId::kNcuts));
   h3v0EtaPtChi2   = addHist(new TH3F("h3v0EtaPtChi2", "#chi^{2};#eta;p_{T} (GeV/#it{c});#chi^{2}", etaNbins, etaMin,
                                      etaMax, ptNbins, ptMin, ptMax, 100, 0., 10.));
   if (isMC) h3v0EtaPtChi2True = addHistClone(h3v0EtaPtChi2, "True");
   h3v0EtaPtRconv = addHist(new TH3F("h3v0EtaPtRconv", "R^{conv};#eta;p_{T} (GeV/#it{c});R^{conv} (cm)", etaNbins,
                                     etaMin, etaMax, ptNbins, ptMin, ptMax, 200, 0., 200.));
   if (isMC) h3v0EtaPtRconvTrue = addHistClone(h3v0EtaPtRconv, "True");
   h3v0EtaPtCPA = addHist(new TH3F("h3v0EtaPtCPA", "CPA distribution;#eta;p_{T} (GeV/#it{c});CPA", etaNbins, etaMin,
                                   etaMax, ptNbins, ptMin, ptMax, 100, -1, 1));
   if (isMC) h3v0EtaPtCPATrue = addHistClone(h3v0EtaPtCPA, "True");
   h3v0EtaPtDist = addHist(new TH3F("h3v0EtaPtDist", "track DCA;#eta;p_{T} (GeV/#it{c});DCA (cm)", etaNbins, etaMin,
                                    etaMax, ptNbins, ptMin, ptMax, 100, 0., 10.));
   if (isMC) h3v0EtaPtDistTrue = addHistClone(h3v0EtaPtDist, "True");
   h3v0EtaPtMassEE = addHist(new TH3F("h3v0EtaPtMassEE", "m_{ee};#eta;p_{T} (GeV/#it{c});m_{ee} (GeV/#it{c}^{2})",
                                      etaNbins, etaMin, etaMax, ptNbins, ptMin, ptMax, 100, 0., 0.3));
   if (isMC) h3v0EtaPtMassEETrue = addHistClone(h3v0EtaPtMassEE, "True");
   h2v0ArmPo = addHist(new TH2F("h2v0ArmPo", "Armenteros plot;#alpha;q*t", 100, -1, 1, 100, 0, 0.3));
   if (isMC) h2v0ArmPoTrue = addHistClone(h2v0ArmPo, "True");
   h3v0EtaPtAsym = addHist(new TH3F("h3v0EtaPtAsym", "Asymmetry;#eta;p_{T} (GeV/#it{c});Asymmetry", etaNbins, etaMin,
                                    etaMax, ptNbins, ptMin, ptMax, 200, 0, 1));
   if (isMC) h3v0EtaPtAsymTrue = addHistClone(h3v0EtaPtAsym, "True");
   h3v0EtaPtCosPsi = addHist(new TH3F("h3v0EtaPtCosPsi", "cos(#psi);#eta;p_{T} (GeV/#it{c});cos(#psi)", etaNbins,
                                      etaMin, etaMax, ptNbins, ptMin, ptMax, 100, -1., 1.));
   if (isMC) h3v0EtaPtCosPsiTrue = addHistClone(h3v0EtaPtCosPsi, "True");
   h2v0EtaPt = addHist(
      new TH2F("h2v0EtaPt", "Conv Sp;#eta;p_{T} (GeV/#it{c})", etaNbins, etaMin, etaMax, ptNbins, ptMin, ptMax));
   if (isMC) h2v0EtaPtTrue = addHistClone(h2v0EtaPt, "True");
   h3v0ConvMap = addHist(new TH3F("h3v0ConvMap", "Conversion map (x,y,z);x (cm);y (cm);z (cm)", 100, -20, 20., 100, -20,
                                  20, 100, -200., 200.));
   if (isMC) h3v0ConvMapTrue = addHistClone(h3v0ConvMap, "True");

   // Cluster selection
   int   eNbins = 100;
   float eMin = 0, eMax = 2.5;
   h3cluEtaECutEff = addHist(new TH3F("h3cluEtaECutEff", "Cluster cut efficiency;#eta;Cut ID", etaNbins, etaMin, etaMax,
                                      eNbins, eMin, eMax, clustCutId::kNcuts, 0., clustCutId::kNcuts));
   h3cluEtaEMult = addHist(new TH3F("h3cluEtaEMult", "Cluster multiplicity;#eta;e (GeV);Mult", etaNbins, etaMin, etaMax,
                                    eNbins, eMin, eMax, 100, 0, 100));
   h3cluEtaETofSigma = addHist(new TH3F("h3cluEtaETofSigma", "Cluster #sigma_{TOF};#eta;e (GeV);#sigma_{TOF}", etaNbins,
                                        etaMin, etaMax, eNbins, eMin, eMax, 200, -10, 10));
   h3cluEtaEDistCPV  = addHist(new TH3F("h3cluEtaEDistCPV", "Cluster DistCPV;#eta;e (GeV);DistCPV (cm)", etaNbins,
                                        etaMin, etaMax, eNbins, eMin, eMax, 100, 0, 30));
   // Inv mass histos
   const int   yNbins = 15, ptNbinsWide = 50;
   const int   mInvNbins = 200;
   const float yMin = -1.5, yMax = 1.5, mInvMax = 1.;
   hMixBinOccupancy  = addHist(new TH1F("hMixBinOccupancy", "hMixBinOccupancy;mixBin;nEvents", nMixBins, 0, nMixBins));
   h3MixBinOccupancy = addHist(new TH3F("h3MixBinOccupancy", "h3MixBinOccupancy;centrality;Zvtx;psiEP", nMixEventCent,
                                        0, nMixEventCent, nMixEventZ, 0, nMixEventZ, nMixEventEP, 0, nMixEventEP));
   h3BinOccupancyRealCalo =
      addHist(new TH3F("h3BinOccupancyRealCalo", "h3BinOccupancyRealCalo;centrality;Zvtx;psiEP", nMixEventCent, 0,
                       nMixEventCent, nMixEventZ, 0, nMixEventZ, nMixEventEP, 0, nMixEventEP));
   h3BinOccupancyMixedCalo =
      addHist(new TH3F("h3BinOccupancyMixedCalo", "h3BinOccupancyMixedCalo;centrality;Zvtx;psiEP", nMixEventCent, 0,
                       nMixEventCent, nMixEventZ, 0, nMixEventZ, nMixEventEP, 0, nMixEventEP));

   // EP resolution
   pR1cent     = addHist(new TProfile("pR1cent", "pR1cent", mCentBins.size() - 1, mCentBins.data()));
   pR1centTrue = addHistClone(pR1cent, "True");

   for (int cen = 0; cen < mHistoCentBins; cen++) {
      std::string centrality = Form("%.0f-%.0f", mAxCent->GetBinLowEdge(cen + 1), mAxCent->GetBinUpEdge(cen + 1));
      h3yPtMinvRealCalo.push_back(addHist(new TH3F(
         Form("h3yPtMinvRealCalo_cen_%s", centrality.c_str()),
         Form("Real inv mass, calorimeter, %s;#it{y};p_{T} (GeV/#it{c});M_{inv} (GeV/#it{c}))", centrality.c_str()),
         yNbins, yMin, yMax, ptNbinsWide, ptMin, ptMax, mInvNbins, 0., mInvMax)));
      if (isMC) {
         h3yPtMinvRealCaloTrueV0.push_back(addHistClone(h3yPtMinvRealCalo.back(), "_TrueV0"));
         h3yPtMinvRealCaloTruePi.push_back(addHistClone(h3yPtMinvRealCalo.back(), "_TruePi"));
         h3yPtMinvRealCaloTruePh.push_back(addHistClone(h3yPtMinvRealCalo.back(), "_TruePh"));
      }
      h3yPtMinvMixedCalo.push_back(addHist(new TH3F(
         Form("h3yPtMinvMixedCalo_cen_%s", centrality.c_str()),
         Form("Mixed inv mass, calorimeter, %s;#it{y};p_{T} (GeV/#it{c});M_{inv} (GeV/#it{c}))", centrality.c_str()),
         yNbins, yMin, yMax, ptNbinsWide, ptMin, ptMax, mInvNbins, 0., mInvMax)));
      h3yPtMinvRealHybrid.push_back(addHist(
         new TH3F(Form("h3yPtMinvRealHybrid_cen_%s", centrality.c_str()),
                  Form("Real inv mass, hybrid, %s;#it{y};p_{T} (GeV/#it{c});M_{inv} (GeV/#it{c}))", centrality.c_str()),
                  yNbins, yMin, yMax, ptNbinsWide, ptMin, ptMax, mInvNbins, 0., mInvMax)));
      if (isMC) {
         h3yPtMinvRealHybridTrueV0.push_back(addHistClone(h3yPtMinvRealHybrid.back(), "_TrueV0"));
         h3yPtMinvRealHybridTruePi.push_back(addHistClone(h3yPtMinvRealHybrid.back(), "_TruePi"));
         h3yPtMinvRealHybridTruePh.push_back(addHistClone(h3yPtMinvRealHybrid.back(), "_TruePh"));
      }
      h3yPtMinvMixedHybrid.push_back(addHist(
         new TH3F(Form("h3yPtMinvMixedHybrid_cen_%s", centrality.c_str()),
                  Form("Mixed inv mass, hybrid, %s;#it{y};p_{T} (GeV/#it{c});M_{inv} (GeV/#it{c})", centrality.c_str()),
                  yNbins, yMin, yMax, ptNbinsWide, ptMin, ptMax, mInvNbins, 0., mInvMax)));
      h3yPtMinvRealConv.push_back(addHist(new TH3F(
         Form("h3yPtMinvRealConv_cen_%s", centrality.c_str()),
         Form("Real inv mass, conversion, %s;#it{y};p_{T} (GeV/#it{c});M_{inv} (GeV/#it{c})", centrality.c_str()),
         yNbins, yMin, yMax, ptNbinsWide, ptMin, ptMax, mInvNbins, 0., mInvMax)));
      if (isMC) {
         h3yPtMinvRealConvTrueV0.push_back(addHistClone(h3yPtMinvRealConv.back(), "_TrueV0"));
         h3yPtMinvRealConvTruePi.push_back(addHistClone(h3yPtMinvRealConv.back(), "_TruePi"));
         h3yPtMinvRealConvTruePh.push_back(addHistClone(h3yPtMinvRealConv.back(), "_TruePh"));
      }
      h3yPtMinvMixedConv.push_back(addHist(new TH3F(
         Form("h3yPtMinvMixedConv_cen_%s", centrality.c_str()),
         Form("Mixed inv mass, conversion, %s;#it{y};p_{T} (GeV/#it{c});M_{inv} (GeV/#it{c})", centrality.c_str()),
         yNbins, yMin, yMax, ptNbinsWide, ptMin, ptMax, mInvNbins, 0., mInvMax)));

      // flow
      p3v1yPtMinvCalo.push_back(addHist(
         new TProfile3D(Form("p3v1yPtMinvCalo_cen_%s", centrality.c_str()),
                        Form("V1, Calo, %s;#it{y};p_{T} (GeV/#it{c});M_{inv} (GeV/#it{c}))", centrality.c_str()),
                        yNbins, yMin, yMax, ptNbinsWide, ptMin, ptMax, mInvNbins, 0., mInvMax)));
      p3v1yPtMinvHybrid.push_back(addHist(
         new TProfile3D(Form("p3v1yPtMinvHybrid_cen_%s", centrality.c_str()),
                        Form("V1, Hybrid, %s;#it{y};p_{T} (GeV/#it{c});M_{inv} (GeV/#it{c}))", centrality.c_str()), 15,
                        -1.5, 1.5, 25, 0, 5, mInvNbins, 0., mInvMax)));
      p3v1yPtMinvConv.push_back(addHist(
         new TProfile3D(Form("p3v1yPtMinvConv_cen_%s", centrality.c_str()),
                        Form("V1, Conv, %s;#it{y};p_{T} (GeV/#it{c});M_{inv} (GeV/#it{c}))", centrality.c_str()),
                        yNbins, yMin, yMax, ptNbinsWide, ptMin, ptMax, mInvNbins, 0., mInvMax)));

      if (isMC) {
         h2yPtPrimPi.push_back(
            addHist(new TH2F(Form("h2yPtPrimPi_cen_%s", centrality.c_str()),
                             Form("Primary Pi0, centrality %s;#it{y};p_{T} (GeV/#it{c})", centrality.c_str()), yNbins,
                             yMin, yMax, ptNbins, ptMin, ptMax)));
         h2yPtPrimPh.push_back(
            addHist(new TH2F(Form("h2yPtPrimPh_cen_%s", centrality.c_str()),
                             Form("Primary Photon, centrality %s;#it{y};p_{T} (GeV/#it{c})", centrality.c_str()),
                             yNbins, yMin, yMax, ptNbins, ptMin, ptMax)));
         p2v1yPtPrimPi.push_back(
            addHist(new TProfile2D(Form("p2v1yPt_cen_%s_Pi", centrality.c_str()),
                                   Form("V1 Pi, True, %s;#it{y};p_{T} (GeV/#it{c}));v_{1}", centrality.c_str()), yNbins,
                                   yMin, yMax, ptNbinsWide, ptMin, ptMax)));
         p2v1yPtPrimPh.push_back(
            addHist(new TProfile2D(Form("p2v1yPt_cen_%s_Ph", centrality.c_str()),
                                   Form("V1 Ph, True, %s;#it{y};p_{T} (GeV/#it{c}));v_{1}", centrality.c_str()), yNbins,
                                   yMin, yMax, ptNbinsWide, ptMin, ptMax)));
         p3v1yPtMinvCaloTruePi.push_back(addHistClone(p3v1yPtMinvCalo.back(), "_TruePi"));
         p3v1yPtMinvHybridTruePi.push_back(addHistClone(p3v1yPtMinvHybrid.back(), "_TruePi"));
         p3v1yPtMinvConvTruePi.push_back(addHistClone(p3v1yPtMinvConv.back(), "_TruePi"));
         p3v1yPtMinvCaloTruePh.push_back(addHistClone(p3v1yPtMinvCalo.back(), "_TruePh"));
         p3v1yPtMinvHybridTruePh.push_back(addHistClone(p3v1yPtMinvHybrid.back(), "_TruePh"));
         p3v1yPtMinvConvTruePh.push_back(addHistClone(p3v1yPtMinvConv.back(), "_TruePh"));
         p3v1yPtMinvCaloRP.push_back(addHistClone(p3v1yPtMinvCalo.back(), "_RP"));
         p3v1yPtMinvHybridRP.push_back(addHistClone(p3v1yPtMinvHybrid.back(), "_RP"));
         p3v1yPtMinvConvRP.push_back(addHistClone(p3v1yPtMinvConv.back(), "_RP"));
         p3v1yPtMinvCaloTruePiRP.push_back(addHistClone(p3v1yPtMinvCalo.back(), "_TruePiRP"));
         p3v1yPtMinvHybridTruePiRP.push_back(addHistClone(p3v1yPtMinvHybrid.back(), "_TruePiRP"));
         p3v1yPtMinvConvTruePiRP.push_back(addHistClone(p3v1yPtMinvConv.back(), "_TruePiRP"));
         p3v1yPtMinvCaloTruePhRP.push_back(addHistClone(p3v1yPtMinvCalo.back(), "_TruePhRP"));
         p3v1yPtMinvHybridTruePhRP.push_back(addHistClone(p3v1yPtMinvHybrid.back(), "_TruePhRP"));
         p3v1yPtMinvConvTruePhRP.push_back(addHistClone(p3v1yPtMinvConv.back(), "_TruePhRP"));
         p2v1yPtPrimPiRP.push_back(addHistClone(p2v1yPtPrimPi.back(), "RP"));
         p2v1yPtPrimPhRP.push_back(addHistClone(p2v1yPtPrimPh.back(), "RP"));
      }
   }

   cout << "MpdConvPi0::UserInit(): complete\n";
}
//--------------------------------------
void MpdConvPi0::ProcessEvent(MpdAnalysisEvent &event)
{
   // Centrality
   mCentrality = event.getCentrTPC();
   mCenBin     = mAxCent->FindBin(mCentrality) - 1;

   if (!selectEvent(event)) return;

   if (mCenBin < 0) mCenBin = 0;
   if (mCenBin >= nMixEventCent) mCenBin = nMixEventCent - 1;
   if (mCenBin >= mHistoCentBins) mCenBin = mHistoCentBins - 1;

   mMpdGlobalTracks = event.fMPDEvent->GetGlobalTracks();
   mKalmanTracks    = event.fTPCKalmanTrack;
   mEMCClusters     = event.fEMCCluster;
   mMCTracks        = event.fMCTrack;
   psiRP            = event.fMCEventHeader->GetRotZ();

   // vtxZ
   mZvtxBin = 0.5 * nMixEventZ * (1 + mPrimaryVertex.Z() / mParams.mZvtxCut);
   if (mZvtxBin < 0) mZvtxBin = 0;
   if (mZvtxBin >= nMixEventZ) mZvtxBin = nMixEventZ - 1;

   // EP
   psiEP  = event.fMpdEP.GetPhiEP_FHCal_F_all();
   psiEPn = event.fMpdEP.GetPhiEP_FHCal_N_all();
   psiEPs = event.fMpdEP.GetPhiEP_FHCal_S_all();
   mEPBin = 0.5 * nMixEventEP * (1 + psiEP / TMath::Pi());
   if (mEPBin < 0) mEPBin = 0;
   if (mEPBin >= nMixEventEP) mEPBin = nMixEventEP - 1;

   pR1cent->Fill(mCentrality, cos(psiEPn - psiEPs));
   pR1centTrue->Fill(mCentrality, cos(psiEP - psiRP));

   int nTracks = mMpdGlobalTracks->GetEntriesFast();
   isGoodTrack.resize(nTracks, true);
   for (int i = 0; i < nTracks; i++) {
      auto *track = static_cast<MpdTrack *>(mMpdGlobalTracks->At(i));
      if (!selectTrack(track)) isGoodTrack.at(i) = false;
   }

   selectConversion(event);

   selectClusters(event);

   processHistograms(event);

   if (isMC) fillMC(event);
}

void MpdConvPi0::fillMC(MpdAnalysisEvent &event)
{
   for (int i = 0; i < mMCTracks->GetEntriesFast(); i++) {
      auto *p = static_cast<MpdMCTrack *>(mMCTracks->At(i));
      if (p->GetStartX() * p->GetStartX() + p->GetStartY() * p->GetStartY() <
          1.) { // make configurable??? what about StartZ?
         TVector3 mom3;
         p->GetMomentum(mom3);
         TLorentzVector mom4;
         float          pt = mom3.Pt(), v1 = cos(mom3.Phi() - psiEP), v1rp = cos(mom3.Phi() - psiRP);
         if (p->GetPdgCode() == 111) {
            float m = TDatabasePDG::Instance()->GetParticle(p->GetPdgCode())->Mass();
            mom4.SetXYZM(mom3.Px(), mom3.Py(), mom3.Pz(), m);
            float y = mom4.Rapidity();
            h2yPtPrimPi.at(mCenBin)->Fill(y, pt);
            p2v1yPtPrimPi.at(mCenBin)->Fill(y, pt, v1);
            p2v1yPtPrimPiRP.at(mCenBin)->Fill(y, pt, v1rp);
         }
         if (p->GetPdgCode() == 22) {
            float m = 0;
            mom4.SetXYZM(mom3.Px(), mom3.Py(), mom3.Pz(), m);
            float y = mom4.Rapidity();
            h2yPtPrimPh.at(mCenBin)->Fill(y, pt);
            p2v1yPtPrimPh.at(mCenBin)->Fill(y, pt, v1);
            p2v1yPtPrimPhRP.at(mCenBin)->Fill(y, pt, v1rp);
         }
      }
   }
}

void MpdConvPi0::Finish()
{
   // Post-scan processing not needed
}

//--------------------------------------
bool MpdConvPi0::selectEvent(MpdAnalysisEvent &event)
{

   hEventCutEff->Fill(eventCutId::kNone);
   // first test if event filled?
   if (!event.fVertex) return false;
   hEventCutEff->Fill(eventCutId::kVertexPresent);

   // Vertex z coordinate
   MpdVertex *vertex = (MpdVertex *)event.fVertex->First();
   vertex->Position(mPrimaryVertex);
   if (applySelection && fabs(mPrimaryVertex.Z()) > mParams.mZvtxCut) return false;
   hVertex->Fill(mPrimaryVertex.Z());
   hEventCutEff->Fill(eventCutId::kVertexZ);

   if (applySelection && (mCenBin < 0 || mCenBin >= mHistoCentBins)) return false;
   hCentrality->Fill(mCentrality);
   hEventCutEff->Fill(eventCutId::kGoodCentrality);

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
      TLorentzVector   mom;
      clu->GetMomentum(mom, &mPrimaryVertex);
      float eta = mom.Eta();
      float e   = Nonlinearity(clu->GetE());

      h3cluEtaECutEff->Fill(eta, e, clustCutId::kNone);

      if (e < mParams.mCluEmin) {
         continue;
      }
      h3cluEtaECutEff->Fill(eta, e, clustCutId::kE);

      int cluMult = clu->GetMultiplicity();
      h3cluEtaEMult->Fill(eta, e, cluMult);
      if (cluMult < mParams.mCluMult) {
         continue;
      }
      h3cluEtaECutEff->Fill(eta, e, clustCutId::kMult);

      double dx   = clu->GetX() - mPrimaryVertex.X();
      double dy   = clu->GetY() - mPrimaryVertex.Y();
      double dz   = clu->GetZ() - mPrimaryVertex.Z();
      double r    = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
      double time = clu->GetTime() - r / 29979245800. * 1.e+9; // Time in ns
      if (isMC) {
         time = smearTime(time, e);
      }

      float cluTofSigma = tofCut(time, e);
      h3cluEtaETofSigma->Fill(eta, e, cluTofSigma);
      if (cluTofSigma < mParams.mCluTofMin || cluTofSigma > mParams.mCluTofMax) {
         continue;
      }

      h3cluEtaECutEff->Fill(eta, e, clustCutId::kTof);

      float l1, l2; // Dispersion axis
      clu->GetLambdas(l1, l2);
      // correct for z
      float izMax = TMath::Abs(32. * clu->GetZ());
      l1          = l1 * 1.6 / (1.6 + 0.0002 * izMax * izMax);
      l2          = l2 * 3.2 / (3.2 + 0.000023 * izMax * izMax * izMax);

      //    if(e>mParams.mCluDispeMin &&  lambdaCut(l1,l2,e) > mParams.mCluDisp){
      //      continue ;
      //    }
      h3cluEtaECutEff->Fill(eta, e, clustCutId::kDisp);

      float pp0      = par0[0] + par0[1] * mPrimaryVertex.Z();
      float pp1      = par1[0] + par1[1] * mPrimaryVertex.Z();
      float z_shift1 = pp0 + pp1 * log(e);

      float cluDistCPV = distCPV(clu->GetDPhi(), clu->GetDZ() + z_shift1, e);
      h3cluEtaEDistCPV->Fill(eta, e, cluDistCPV, e);
      if (cluDistCPV < mParams.mCluCPV) {
         continue;
      }
      h3cluEtaECutEff->Fill(eta, e, clustCutId::kCPV);

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
         float          y = sum.Rapidity(), pt = sum.Pt(), m = sum.M(), v1 = cos(sum.Phi() - psiEP),
               v1rp = cos(sum.Phi() - psiRP);
         h3BinOccupancyRealCalo->Fill(mCenBin, mZvtxBin, mEPBin);
         h3yPtMinvRealCalo.at(mCenBin)->Fill(y, pt, m);
         p3v1yPtMinvCalo.at(mCenBin)->Fill(y, pt, m, v1);
         p3v1yPtMinvCaloRP.at(mCenBin)->Fill(y, pt, m, v1rp);
         long int ip = MpdV0Maker::FindCommonParent(mClusters[i].primary(), mClusters[j].primary(), mMCTracks);
         if (ip >= 0) {
            int pdg = abs(static_cast<MpdMCTrack *>(mMCTracks->At(ip))->GetPdgCode());
            h3yPtMinvRealCaloTrueV0.at(mCenBin)->Fill(y, pt, m);
            if (pdg == 111) {
               h3yPtMinvRealCaloTruePi.at(mCenBin)->Fill(y, pt, m);
               p3v1yPtMinvCaloTruePi.at(mCenBin)->Fill(y, pt, m, v1);
               p3v1yPtMinvCaloTruePiRP.at(mCenBin)->Fill(y, pt, m, v1rp);
            }
            if (pdg == 22) {
               h3yPtMinvRealCaloTruePh.at(mCenBin)->Fill(y, pt, m);
               p3v1yPtMinvCaloTruePh.at(mCenBin)->Fill(y, pt, m, v1);
               p3v1yPtMinvCaloTruePhRP.at(mCenBin)->Fill(y, pt, m, v1rp);
            }
         }
      }
   }

   for (int i = 0; i < nClu; i++) {
      for (int j = 0; j < nV0; j++) {
         if (!TestHybrid(mClusters[i], mV0[j])) continue;
         TLorentzVector sum = mClusters[i] + mV0[j];
         float          y = sum.Rapidity(), pt = sum.Pt(), m = sum.M(), v1 = cos(sum.Phi() - psiEP),
               v1rp = cos(sum.Phi() - psiRP);
         h3yPtMinvRealHybrid.at(mCenBin)->Fill(y, pt, m);
         p3v1yPtMinvHybrid.at(mCenBin)->Fill(y, pt, m, v1);
         p3v1yPtMinvHybridRP.at(mCenBin)->Fill(y, pt, m, v1rp);
         long int ip = MpdV0Maker::FindCommonParent(mClusters[i].primary(), mV0[j].primary(), mMCTracks);
         if (ip >= 0) {
            int pdg = abs(static_cast<MpdMCTrack *>(mMCTracks->At(ip))->GetPdgCode());
            h3yPtMinvRealHybridTrueV0.at(mCenBin)->Fill(y, pt, m);
            if (pdg == 111) {
               h3yPtMinvRealHybridTruePi.at(mCenBin)->Fill(y, pt, m);
               p3v1yPtMinvHybridTruePi.at(mCenBin)->Fill(y, pt, m, v1);
               p3v1yPtMinvHybridTruePiRP.at(mCenBin)->Fill(y, pt, m, v1rp);
            }
            if (pdg == 22) {
               h3yPtMinvRealHybridTruePh.at(mCenBin)->Fill(y, pt, m);
               p3v1yPtMinvHybridTruePh.at(mCenBin)->Fill(y, pt, m, v1);
               p3v1yPtMinvHybridTruePhRP.at(mCenBin)->Fill(y, pt, m, v1rp);
            }
         }
      }
   }

   for (int i = 0; i < nV0 - 1; i++) {
      for (int j = i + 1; j < nV0; j++) {
         TLorentzVector sum = mV0[i] + mV0[j];
         float          y = sum.Rapidity(), pt = sum.Pt(), m = sum.M(), v1 = cos(sum.Phi() - psiEP),
               v1rp = cos(sum.Phi() - psiRP);
         h3yPtMinvRealConv.at(mCenBin)->Fill(y, pt, m);
         p3v1yPtMinvConv.at(mCenBin)->Fill(y, pt, m, v1);
         p3v1yPtMinvConvRP.at(mCenBin)->Fill(y, pt, m, v1rp);
         long int ip = MpdV0Maker::FindCommonParent(mV0[i].primary(), mV0[j].primary(), mMCTracks);
         if (ip >= 0) {
            int pdg = abs(static_cast<MpdMCTrack *>(mMCTracks->At(ip))->GetPdgCode());
            h3yPtMinvRealConvTrueV0.at(mCenBin)->Fill(y, pt, m);
            if (pdg == 111) {
               h3yPtMinvRealConvTruePi.at(mCenBin)->Fill(y, pt, m);
               p3v1yPtMinvConvTruePi.at(mCenBin)->Fill(y, pt, m, v1);
               p3v1yPtMinvConvTruePiRP.at(mCenBin)->Fill(y, pt, m, v1rp);
            }
            if (pdg == 22) {
               h3yPtMinvRealConvTruePh.at(mCenBin)->Fill(y, pt, m);
               p3v1yPtMinvConvTruePh.at(mCenBin)->Fill(y, pt, m, v1);
               p3v1yPtMinvConvTruePhRP.at(mCenBin)->Fill(y, pt, m, v1rp);
            }
         }
      }
   }

   // Mixed
   // calculate bin from zVertex-centrality-reaction plane
   int mixBin = mCenBin * nMixEventEP * nMixEventZ + mEPBin * nMixEventZ + mZvtxBin;
   hMixBinOccupancy->Fill(mixBin);
   h3MixBinOccupancy->Fill(mCenBin, mZvtxBin, mEPBin);

   for (auto &mixEventClusters : mMixEventsClusters.at(mixBin)) {
      for (auto &vm : mixEventClusters) {
         for (auto &v : mClusters) {
            TLorentzVector sum = v + vm;
            h3BinOccupancyMixedCalo->Fill(mCenBin, mZvtxBin, mEPBin);
            h3yPtMinvMixedCalo.at(mCenBin)->Fill(sum.Rapidity(), sum.Pt(), sum.M());
         }
      }
      for (auto &vm : mixEventClusters) {
         for (auto &v : mV0) {
            TLorentzVector sum = v + vm;
            h3yPtMinvMixedHybrid.at(mCenBin)->Fill(sum.Rapidity(), sum.Pt(), sum.M());
         }
      }
   }
   for (auto &mixEventV0 : mMixEventsV0.at(mixBin)) {
      for (auto &vm : mixEventV0) {
         for (auto &v : mClusters) {
            TLorentzVector sum = v + vm;
            h3yPtMinvMixedHybrid.at(mCenBin)->Fill(sum.Rapidity(), sum.Pt(), sum.M());
         }
      }
      for (auto &vm : mixEventV0) {
         for (auto &v : mV0) {
            TLorentzVector sum = v + vm;
            h3yPtMinvMixedConv.at(mCenBin)->Fill(sum.Rapidity(), sum.Pt(), sum.M());
         }
      }
   }

   // Append new particles to queue and remove those at the beginning
   if (mClusters.size() > 0) mMixEventsClusters.at(mixBin).emplace_back(mClusters);
   if (mMixEventsClusters.at(mixBin).size() > mMaxMixSize) mMixEventsClusters.at(mixBin).pop_front();
   if (mV0.size() > 0) mMixEventsV0.at(mixBin).emplace_back(mV0);
   if (mMixEventsV0.at(mixBin).size() > mMaxMixSize) mMixEventsV0.at(mixBin).pop_front();
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
   h3trackEtaPtCutEff->Fill(eta, pt, trackCutId::kNone);

   int nHits = mpdtrack->GetNofHits();
   if (isElectron) h3trackEtaPtNhitsTrue->Fill(eta, pt, nHits);
   if (applySelection && nHits < mParams.mNofHitsCut) return false;
   h3trackEtaPtNhits->Fill(eta, pt, nHits);
   h3trackEtaPtCutEff->Fill(eta, pt, trackCutId::kNhits);

   if (isElectron) h2trackEtaPtTrue->Fill(eta, pt);
   if (applySelection && !(fabs(eta) < mParams.mEtaCut && pt > mParams.mPtminCut)) return false;
   h2trackEtaPt->Fill(eta, pt);
   h3trackEtaPtCutEff->Fill(eta, pt, trackCutId::kEtaPt);

   float probEl = mpdtrack->GetPidProbElectron(); // TPC only or combined??? (actually none works)
   if (isElectron) h3trackEtaPtProbElTrue->Fill(eta, pt, probEl);
   if (applySelection && probEl < mParams.mProbElCut) return false;
   h3trackEtaPtProbEl->Fill(eta, pt, probEl);
   h3trackEtaPtCutEff->Fill(eta, pt, trackCutId::kProbElectron);

   bool  hasTofHit   = (fabs(mpdtrack->GetTofDphiSigma()) < 3.0 && fabs(mpdtrack->GetTofDzSigma()) < 3.0);
   bool  tofElectron = false, dEdxElectron = false;
   float nSigma_dEdx = mpdtrack->GetTPCNSigma(kEl);
   float nSigma_beta = mpdtrack->GetTOFNSigma(kEl);
   if (isElectron) { // same for true electron tracks
      h3trackEtaPtNsigDEdxTrue->Fill(eta, pt, nSigma_dEdx);
      h3trackEtaPtNsigBetaTrue->Fill(eta, pt, nSigma_beta);
      h2trackNsigDEdxNsigBetaTrue->Fill(nSigma_dEdx, nSigma_beta);
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
   h3trackEtaPtNsigDEdx->Fill(eta, pt, nSigma_dEdx);
   h3trackEtaPtNsigBeta->Fill(eta, pt, nSigma_beta);
   h2trackNsigDEdxNsigBeta->Fill(nSigma_dEdx, nSigma_beta);
   h3trackEtaPtCutEff->Fill(eta, pt, trackCutId::kTpcTofPid);

   return true;
}

bool MpdConvPi0::selectV0(MpdV0 *v0)
{
   // Check V0
   float pt = v0->Pt(), eta = v0->Eta();
   h3v0EtaPtCutEff->Fill(eta, pt, V0CutId::kNone);

   if (applySelection && mParams.mUseBDT) {
      if (v0->getBDTValue() > mParams.mBDTCut) {
         h3v0EtaPtCutEff->Fill(eta, pt, V0CutId::kBDT);
         if (mParams.mUseBDTRegP) {
            double scale = v0->getBDTMom() / v0->P();
            (*v0) *= scale;
         }
         return true;
      } else
         return false;
   }

   if (!(isGoodTrack.at(v0->getTr1()) && isGoodTrack.at(v0->getTr2()))) return false;
   h3v0EtaPtCutEff->Fill(eta, pt, V0CutId::kGoodTracks);

   if (applySelection && pt < 0.005) return false; // to avoid fpe
   h3v0EtaPtCutEff->Fill(eta, pt, V0CutId::kNonZeroPt);

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
      h3v0EtaPtChi2True->Fill(eta, pt, chi2);
   }
   if (applySelection && chi2 > mParams.mChi2Cut) return false;
   h3v0EtaPtChi2->Fill(eta, pt, chi2);
   h3v0EtaPtCutEff->Fill(eta, pt, V0CutId::kChi2);

   float x0, y0, z0;
   v0->getVertex(x0, y0, z0);

   h3v0ConvMap->Fill(x0, y0, z0);
   if (isTrue) {
      h3v0ConvMapTrue->Fill(x0, y0, z0);
   }

   float rConv = v0->getRconv();
   if (isTrue) { // same for true electrontracks
      h3v0EtaPtRconvTrue->Fill(eta, pt, rConv);
   }
   if (applySelection && (rConv < mParams.mMinR2Cut || rConv > mParams.mMaxR2Cut)) return false;
   h3v0EtaPtRconv->Fill(eta, pt, rConv);
   h3v0EtaPtCutEff->Fill(eta, pt, V0CutId::kRconv);

   float cpa = v0->getCPA();
   if (isTrue) { // same for true electrontracks
      h3v0EtaPtCPATrue->Fill(eta, pt, cpa);
   }
   if (applySelection && cpa > mParams.mCPACut) return false;
   h3v0EtaPtCPA->Fill(eta, pt, cpa);
   h3v0EtaPtCutEff->Fill(eta, pt, V0CutId::kCPA);

   float dist = v0->getDaughterDCA();

   if (isTrue) { // same for true electrontracks
      h3v0EtaPtDistTrue->Fill(eta, pt, dist);
   }
   if (applySelection && dist > mParams.mDistCut) return false;
   h3v0EtaPtDist->Fill(eta, pt, dist);
   h3v0EtaPtCutEff->Fill(eta, pt, V0CutId::kDist);

   float mass = v0->getMass();
   if (isTrue) { // same for true electrontracks
      h3v0EtaPtMassEETrue->Fill(eta, pt, mass);
   }
   if (applySelection && mass > mParams.mMassCut) return false;
   h3v0EtaPtMassEE->Fill(eta, pt, mass);
   h3v0EtaPtCutEff->Fill(eta, pt, V0CutId::kMass);

   float alpha, qt;
   v0->getArmenteros(alpha, qt);
   if (isTrue) {
      h2v0ArmPoTrue->Fill(alpha, qt);
   }
   if (applySelection && !ArmenterosQtCut(v0)) {
      // return false;
   }
   h2v0ArmPo->Fill(alpha, qt);
   h3v0EtaPtCutEff->Fill(eta, pt, V0CutId::kQt);

   // Asymmetry cut
   float asym1, asym2;
   v0->getAsymmetry(asym1, asym2);
   if (isTrue) {
      h3v0EtaPtAsymTrue->Fill(eta, pt, asym1);
      h3v0EtaPtAsymTrue->Fill(eta, pt, asym2);
   }
   if (applySelection && !(AsymmetryCut(asym1, pt) && AsymmetryCut(asym2, pt))) {
      return false;
   }
   h3v0EtaPtAsym->Fill(eta, pt, asym1, pt);
   h3v0EtaPtAsym->Fill(eta, pt, asym2, pt);
   h3v0EtaPtCutEff->Fill(eta, pt, V0CutId::kAsymmetry);

   float cospsi = v0->getCospsi();
   if (isTrue) { // same for true electrontracks
      h3v0EtaPtCosPsiTrue->Fill(eta, pt, cospsi);
   }
   if (applySelection && 1 - fabs(cospsi) > mParams.mCosPsiCut) return false;
   h3v0EtaPtCosPsi->Fill(eta, pt, cospsi, pt);
   h3v0EtaPtCutEff->Fill(eta, pt, V0CutId::kCosPsi);

   h3v0EtaPtCutEff->Fill(eta, pt, V0CutId::kPhotonAOD);

   h2v0EtaPt->Fill(eta, pt);
   if (isTrue) {
      h2v0EtaPtTrue->Fill(eta, pt);
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

float MpdConvPi0::smearTime(float time, float E) const
{
   float sigma = 1.86166 * TMath::Exp(-E / 0.0259728) + 0.347552; // resolution in ns
   return gRandom->Gaus(time, sigma);
}

float MpdConvPi0::Nonlinearity(float oldE) const
{

   float x = TMath::Min(oldE, 2.5f);
   return 2.9411765 * oldE / (0.97630219 + 7.194380e-002 * x - 4.491255e-002 * x * x + 8.362250e-003 * x * x * x);
}
