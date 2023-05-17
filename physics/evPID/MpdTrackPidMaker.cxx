#include <iostream>
#include "MpdTrackPidMaker.h"
#include "MpdTrack.h"
#include "MpdVertex.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdKalmanHit.h"
#include "MpdKalmanFilter.h"
#include "MpdTofMatching.h"
#include "MpdTofMatchingData.h"
#include "MpdEmcClusterKI.h"
//#include "TRandom.h"
#include "TFile.h"

ClassImp(MpdTrackPidMaker);

using namespace std;
MpdTrackPidMaker::MpdTrackPidMaker(const char *name, const char *outputName) : MpdAnalysisTask(name, outputName)
{
   // Create PID engine
}

void MpdTrackPidMaker::UserInit()
{
   cout << "[MpdTrackPidMaker]: Initialization ... " << endl;

   // mParams.ReadFromFile(mParamConfig);
   // mParams.Print();

   // Prepare histograms etc.
   fOutputList = new TList();
   fOutputList->SetOwner(kTRUE);

   TH1::AddDirectory(kFALSE); // sets a global switch disabling the reference to histos in gROOT and their overwriting

   // General QA
   mhEvents = new TH1F("hEvents", "Number of events", 10, 0., 10.);
   fOutputList->Add(mhEvents);
   mhVertex = new TH1F("hVertex", "Event vertex distribution", 100, -200., 200.);
   fOutputList->Add(mhVertex);
   mhCentrality = new TH1F("hCentrality", "Centrality distribution", 100, 0., 100.);
   fOutputList->Add(mhCentrality);
   mhDCAXs = new TH2F("DCAXs", "n-sigma DCAX", 100, 0., 3., 100, -10., 10.);
   fOutputList->Add(mhDCAXs);
   mhDCAYs = new TH2F("DCAYs", "n-sigma DCAY", 100, 0., 3., 100, -10., 10.);
   fOutputList->Add(mhDCAYs);
   mhDCAZs = new TH2F("DCAZs", "n-sigma DCAZ", 100, 0., 3., 100, -10., 10.);
   fOutputList->Add(mhDCAZs);
   mhTPCEl = new TH2F("TPCEl", "n-sigma TPC e-ID", 100, 0., 3., 100, -10., 10.);
   fOutputList->Add(mhTPCEl);
   mhTPCPi = new TH2F("TPCPi", "n-sigma TPC pi-ID", 100, 0., 3., 100, -10., 10.);
   fOutputList->Add(mhTPCPi);
   mhTPCK = new TH2F("TPCK", "n-sigma TPC K-ID", 100, 0., 3., 100, -10., 10.);
   fOutputList->Add(mhTPCK);
   mhTPCP = new TH2F("TPCP", "n-sigma TPC p-ID", 100, 0., 3., 100, -10., 10.);
   fOutputList->Add(mhTPCP);
   mhTOFsdPhi = new TH2F("TOFsdPhi", "n-sigma dPhi matching to TOF", 100, 0., 3., 100, -10., 10.);
   fOutputList->Add(mhTOFsdPhi);
   mhTOFsdZed = new TH2F("TOFsdZed", "n-sigma dZed matching to TOF", 100, 0., 3., 100, -10., 10.);
   fOutputList->Add(mhTOFsdZed);
   mhTOFEl = new TH2F("TOFEl", "n-sigma TOF e-ID", 100, 0., 3., 100, -10., 10.);
   fOutputList->Add(mhTOFEl);
   mhTOFPi = new TH2F("TOFPi", "n-sigma TOF pi-ID", 100, 0., 3., 100, -10., 10.);
   fOutputList->Add(mhTOFPi);
   mhTOFK = new TH2F("TOFK", "n-sigma TOF K-ID", 100, 0., 3., 100, -10., 10.);
   fOutputList->Add(mhTOFK);
   mhTOFP = new TH2F("TOFP", "n-sigma TOF p-ID", 100, 0., 3., 100, -10., 10.);
   fOutputList->Add(mhTOFP);
   mhECALsdPhi = new TH2F("ECALsdPhi", "n-sigma dPhi matching to ECAL", 100, 0., 3., 100, -10., 10.);
   fOutputList->Add(mhECALsdPhi);
   mhECALsdZed = new TH2F("ECALsdZed", "n-sigma dZed matching to ECAL", 100, 0., 3., 100, -10., 10.);
   fOutputList->Add(mhECALsdZed);
   mhEP = new TH2F("EP", "E/p", 100, 0., 3., 300, 0., 3.);
   fOutputList->Add(mhEP);
   mhEM2 = new TH2F("EM2", "ECAL mass2", 100, 0., 3., 300, -0.5, 2.5);
   fOutputList->Add(mhEM2);

   // DCAs
   cout << "[MpdPairKK]: Reading DCA parameterizations ... " << endl;

   // Read-out DCA parameterization
   dcaFile = new TFile("DCAs.root", "READ");

   for (Int_t etab = 0; etab < neta_bins; etab++) {
      for (Int_t centb = 0; centb < ncent_bins; centb++) {
         f_dca_xy[etab][centb] = (TF1 *)dcaFile->Get(Form("dcaxy_fitf_eta_%d_cent_%d", etab, centb));
         f_dca_z[etab][centb]  = (TF1 *)dcaFile->Get(Form("dcaz_fitf_eta_%d_cent_%d", etab, centb));
      }
   }

   TDatime tim;
   RND.SetSeed(tim.GetTime());

   pKF = MpdKalmanFilter::Instance("KF", "KF");

   cout << "[MpdTrackPidMaker]: Initialization done " << endl << endl;
}
//--------------------------------------
void MpdTrackPidMaker::ProcessEvent(MpdAnalysisEvent &event)
{
   // fRunNumber = event.getRunNUmber;
   // if(fParams){ //create parameterizations and readparameters from DB
   // }

   if (!selectEvent(event)) {
      return;
   }

   TClonesArray *mpdGlobalTracks = event.fMPDEvent->GetGlobalTracks();
   TClonesArray *mKalmanTracks   = event.fTPCKalmanTrack;
   TClonesArray *mpdTofMatching  = event.fTOFMatching;
   TObjArray    *mpdEMCClusters  = event.fEMCCluster;

   int nTracks = mpdGlobalTracks->GetEntriesFast();

   for (long int i = 0; i < nTracks; i++) {
      MpdTrack          *mpdtrack = (MpdTrack *)mpdGlobalTracks->UncheckedAt(i);
      MpdTpcKalmanTrack *tr       = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(i);

      // n-sigma DCA
      int cent_bin = -1;
      cent_bin     = int(cen / (cent_max - cent_min) * float(ncent_bins));
      if (cent_bin == 9) cent_bin = 8; // 90-91% -> use parameterization for 80-90%
      if (cent_bin < 0 || cent_bin > ncent_bins - 1) {
         cout << "[MpdTrackPidMaker]: Wrong centrality bin " << cent_bin << " for centrality " << cen
              << " ... skipping event " << endl;
         continue;
      }

      int eta_bin = -1;
      eta_bin     = int((mpdtrack->GetEta() + 1.5) / (eta_max - eta_min) * float(neta_bins));
      if (eta_bin < 0) eta_bin = 0.;
      if (eta_bin > neta_bins - 1)
         eta_bin = neta_bins - 1; // parameterization is valid at |eta| < 1.5 but also used outside

      TF1 sigma_fit_XY = *f_dca_xy[eta_bin][cent_bin];
      TF1 sigma_fit_Z  = *f_dca_z[eta_bin][cent_bin];

      // use lower/upper limits outside of parameterization ranges in pT
      float pt     = fabs(mpdtrack->GetPt());
      float pt_dca = pt;
      if (pt_dca < 0.05) pt_dca = 0.05;
      double xfirst, xlast;
      sigma_fit_Z.GetRange(xfirst, xlast);
      if (pt_dca > xlast) pt_dca = xlast;

      float sigma_exp_xy = sigma_fit_XY(pt_dca);
      float sigma_exp_z  = sigma_fit_Z(pt_dca);

      float dcax_sig, dcay_sig, dcaz_sig;

      if (sigma_exp_xy != 0) {
         dcax_sig = mpdtrack->GetDCAX() / sigma_exp_xy;
         dcay_sig = mpdtrack->GetDCAY() / sigma_exp_xy;
      } else {
         dcax_sig = -999;
         dcay_sig = -999;
      }

      if (sigma_exp_z != 0) {
         dcaz_sig = mpdtrack->GetDCAZ() / sigma_exp_z;
      } else {
         dcaz_sig = -999;
      }

      mpdtrack->SetNSigmaDCAx(dcax_sig);
      mpdtrack->SetNSigmaDCAy(dcay_sig);
      mpdtrack->SetNSigmaDCAz(dcaz_sig);

      mhDCAXs->Fill(pt, mpdtrack->GetNSigmaDCAx());
      mhDCAYs->Fill(pt, mpdtrack->GetNSigmaDCAy());
      mhDCAZs->Fill(pt, mpdtrack->GetNSigmaDCAz());

      // n-sigma dE/dx
      float dEdxsigmas[kNPidLines];
      for (int l = 0; l < 8; l++) {
         dEdxsigmas[l] = -999;
      }

      float pmom = sqrt(mpdtrack->GetPt() * mpdtrack->GetPt() + mpdtrack->GetPz() * mpdtrack->GetPz());

      dEdxsigmas[kEl] = dEdx_sigma_El(tr->GetDedx(), pmom);
      dEdxsigmas[kPi] = dEdx_sigma_Pi(tr->GetDedx(), pmom);
      dEdxsigmas[kK]  = dEdx_sigma_K(tr->GetDedx(), pmom);
      dEdxsigmas[kP]  = dEdx_sigma_P(tr->GetDedx(), pmom);

      mpdtrack->SetTPCNSigma(dEdxsigmas);

      mhTPCEl->Fill(pmom, mpdtrack->GetTPCNSigma(kEl));
      mhTPCPi->Fill(pmom, mpdtrack->GetTPCNSigma(kPi));
      mhTPCK->Fill(pmom, mpdtrack->GetTPCNSigma(kK));
      mhTPCP->Fill(pmom, mpdtrack->GetTPCNSigma(kP));

      // n-sigma matching to TOF
      float tofsdPhi = -999;
      float tofsdZed = -999;

      if (mpdtrack->GetTofFlag() == 2 || mpdtrack->GetTofFlag() == 6) {
         int matchingIndex = -1;

         for (Int_t l = 0; l < mpdTofMatching->GetEntries(); l++) {
            MpdTofMatchingData *Matching = (MpdTofMatchingData *)mpdTofMatching->At(l);
            if (Matching->GetKFTrackIndex() == i) {
               matchingIndex = l;
               break;
            }
         }

         if (matchingIndex > 0) {
            MpdTofMatchingData *Matching = (MpdTofMatchingData *)mpdTofMatching->At(matchingIndex);
            tofsdPhi = TofdPhiMatch(mpdtrack->GetCharge(), fabs(mpdtrack->GetPt()), Matching->GetdPhi());
            tofsdZed = TofdZedMatch(fabs(mpdtrack->GetPt()), Matching->GetdZed());
         }
      } // flag 2||6

      mpdtrack->SetTofMatching(tofsdPhi, tofsdZed);

      mhTOFsdPhi->Fill(fabs(mpdtrack->GetPt()), mpdtrack->GetTofDphiSigma());
      mhTOFsdZed->Fill(fabs(mpdtrack->GetPt()), mpdtrack->GetTofDzSigma());

      // n-sigma TOF-ID
      float tofsigmas[kNPidLines];
      for (int l = 0; l < 8; l++) {
         tofsigmas[l] = -999;
      }

      if (mpdtrack->GetTofFlag() == 2 || mpdtrack->GetTofFlag() == 6) {
         tofsigmas[kEl] = Beta_sigma_El(mpdtrack->GetTofBeta(), pmom);
         tofsigmas[kPi] = Beta_sigma_Pi(mpdtrack->GetTofBeta(), pmom);
         tofsigmas[kK]  = Beta_sigma_K(mpdtrack->GetTofBeta(), pmom);
         tofsigmas[kP]  = Beta_sigma_P(mpdtrack->GetTofBeta(), pmom);
      } // flag 2||6

      mpdtrack->SetTOFNSigma(tofsigmas);

      mhTOFEl->Fill(pmom, mpdtrack->GetTOFNSigma(kEl));
      mhTOFPi->Fill(pmom, mpdtrack->GetTOFNSigma(kPi));
      mhTOFK->Fill(pmom, mpdtrack->GetTOFNSigma(kK));
      mhTOFP->Fill(pmom, mpdtrack->GetTOFNSigma(kP));

      // n-sigma matching ECAL
      float     emcsdPhi = -999;
      float     emcsdZed = -999;
      float     ETr      = -999;
      float     TTr      = -999;
      float     LgTr     = -999;
      int       indTr    = -999;
      float     dDmin, dD;
      float     xCl, yCl, zCl, rCl, phiCl, mass2;
      float     xTr, yTr, zTr, phiTr, emcdPhi, emcdZed, rE;
      int       iok;
      float     rEmin   = 181.7;
      float     rEmax   = 188.1;
      const int rE_bins = 20;
      float     xyzTr[rE_bins][4];

      // slow procedure -> run for tracks identified as electrons in the TPC
      // if ( fabs(mpdtrack -> GetTPCNSigma(kEl)) < 3.0 )
      if (fabs(mpdtrack->GetTPCNSigma(kEl)) < 3.0 && fabs(mpdtrack->GetPt()) > 0.09) {
         if (fabs((*tr->GetParamAtHit())(1, 0)) < 310 && mpdEMCClusters->GetEntries() > 0) {
            MpdKalmanHit hEnd;
            hEnd.SetType(MpdKalmanHit::kFixedR);
            dDmin = 1e9;

            for (int l = 0; l < rE_bins; l++) {
               xyzTr[l][0] = -999;
               xyzTr[l][1] = -999;
               xyzTr[l][2] = -999;
               xyzTr[l][3] = -999;

               rCl = rEmin + (rEmax - rEmin) / float(rE_bins) * (0.5 + l);
               hEnd.SetPos(rCl);

               MpdTpcKalmanTrack tr1(*tr); // new
               tr1.SetParam(*tr1.GetParamAtHit());
               tr1.SetParamNew(*tr1.GetParamAtHit());
               tr1.SetWeight(*tr1.GetWeightAtHit());
               tr1.SetPos(tr1.GetPosAtHit());
               tr1.SetPosNew(tr1.GetPos());
               tr1.SetLength(tr1.GetLengAtHit());

               iok = pKF->PropagateToHit(&tr1, &hEnd, kTRUE);

               if (iok) {
                  phiTr       = tr1.GetParamNew(0) / tr1.GetPosNew();
                  xyzTr[l][0] = tr1.GetPosNew() * cos(phiTr);
                  xyzTr[l][1] = tr1.GetPosNew() * sin(phiTr);
                  xyzTr[l][2] = tr1.GetParamNew(1);
                  xyzTr[l][3] = tr1.GetParamNew(0) / tr1.GetPosNew();
                  // cout<<"Track: "<<i<<" attempt: "<<l<<"  ||  "<<xyzTr[l][0]<<" "<<xyzTr[l][1]<<" "<<xyzTr[l][2]<<"
                  // "<<xyzTr[l][3]<<" R: "<<rCl<<endl;
               }

            } // l

            for (Int_t clu = 0; clu < mpdEMCClusters->GetEntries(); ++clu) {
               MpdEmcClusterKI *EMCCluster = (MpdEmcClusterKI *)mpdEMCClusters->At(clu);

               if (EMCCluster->GetE() / 0.324 < 0.05) continue;
               if (EMCCluster->GetMultiplicity() < 2) continue;
               xCl   = EMCCluster->GetX();
               yCl   = EMCCluster->GetY();
               zCl   = EMCCluster->GetZ();
               rCl   = sqrt(xCl * xCl + yCl * yCl);
               phiCl = TMath::ATan2(yCl, xCl);

               int tr_bin = int((rCl - rEmin) / (rEmax - rEmin) * float(rE_bins));
               if (tr_bin > rE_bins) tr_bin = rE_bins;
               if (tr_bin < 0) tr_bin = 0;

               xTr   = xyzTr[tr_bin][0];
               yTr   = xyzTr[tr_bin][1];
               zTr   = xyzTr[tr_bin][2];
               phiTr = xyzTr[tr_bin][3];

               dD = sqrt(pow(xTr - xCl, 2) + pow(yTr - yCl, 2) + pow(zTr - zCl, 2));

               // cout<<"Track: "<<i<<" Clu: "<<clu<<"  ||  "<<xCl<<" "<<yCl<<" "<<zCl<<" "<<phiCl<<" R: "<<rCl<<"
               // "<<dD<<endl;

               if (dD < dDmin && xTr > -999) {
                  dDmin = dD;
                  indTr = clu;
               } // min
            }    // clu

            if (indTr > 0) {
               // Refit for the best match
               MpdEmcClusterKI *EMCCluster = (MpdEmcClusterKI *)mpdEMCClusters->At(indTr);
               xCl                         = EMCCluster->GetX();
               yCl                         = EMCCluster->GetY();
               zCl                         = EMCCluster->GetZ();
               rCl                         = sqrt(xCl * xCl + yCl * yCl);
               phiCl                       = TMath::ATan2(yCl, xCl);

               hEnd.SetPos(rCl);

               MpdTpcKalmanTrack tr1(*tr); // new
               tr1.SetParam(*tr1.GetParamAtHit());
               tr1.SetParamNew(*tr1.GetParamAtHit());
               tr1.SetWeight(*tr1.GetWeightAtHit());
               tr1.SetPos(tr1.GetPosAtHit());
               tr1.SetPosNew(tr1.GetPos());
               tr1.SetLength(tr1.GetLengAtHit());

               iok = pKF->PropagateToHit(&tr1, &hEnd, kTRUE);

               if (iok) {
                  phiTr = tr1.GetParamNew(0) / tr1.GetPosNew();
                  xTr   = tr1.GetPosNew() * cos(phiTr);
                  yTr   = tr1.GetPosNew() * sin(phiTr);
                  zTr   = tr1.GetParamNew(1);

                  ETr = EMCCluster->GetE();
                  TTr = EMCCluster->GetTime();

                  emcdPhi = phiTr - phiCl;
                  if (emcdPhi < -3.14159265359) emcdPhi = emcdPhi + 2 * 3.14159265359;
                  if (emcdPhi > 3.14159265359) emcdPhi = emcdPhi - 2 * 3.14159265359;

                  emcdZed = zTr - zCl;

                  LgTr = tr1.GetLength();
               } // refit
            }    // indTr > 0
         }       // 310
      }          // el

      emcsdPhi = EmcdPhiMatch(mpdtrack->GetCharge(), fabs(mpdtrack->GetPt()), emcdPhi);
      emcsdZed = EmcdZedMatch(fabs(mpdtrack->GetPt()), emcdZed);
      mpdtrack->SetEcalMatching(emcsdPhi, emcsdZed);

      mhECALsdPhi->Fill(fabs(mpdtrack->GetPt()), mpdtrack->GetECALDphiSigma());
      mhECALsdZed->Fill(fabs(mpdtrack->GetPt()), mpdtrack->GetECALDzSigma());

      if (TTr > -999) TTr = RND.Gaus(TTr, 0.5); // 500ps time resolution
      mpdtrack->SetEt(TTr);

      mpdtrack->SetEp(ETr);
      mpdtrack->SetEl(LgTr);
      mpdtrack->SetEi(indTr);

      if (fabs(mpdtrack->GetECALDphiSigma()) < 3.0 && fabs(mpdtrack->GetECALDzSigma()) < 3.0) {
         if (fabs(mpdtrack->GetTOFNSigma(kEl)) < 3.0) {
            mhEP->Fill(fabs(mpdtrack->GetPt()), (mpdtrack->GetEp()) / 0.324 / pmom);
            mass2 = pmom * pmom *
                    (mpdtrack->GetEt() * mpdtrack->GetEt() * 30 * 30 / mpdtrack->GetEl() / mpdtrack->GetEl() - 1);
            mhEM2->Fill(fabs(mpdtrack->GetPt()), mass2);
         }
      }

   } // ntracks
}

void MpdTrackPidMaker::Finish()
{
   // Post-scan processing not needed
}

//--------------------------------------
bool MpdTrackPidMaker::selectEvent(MpdAnalysisEvent &event)
{
   mhEvents->Fill(0.5);

   if (!event.fVertex) { // if even vertex not filled, skip event
      return false;
   }

   mhEvents->Fill(1.5);

   MpdVertex *vertex = (MpdVertex *)event.fVertex->First();
   vertex->Position(mPrimaryVertex);

   cen = event.getCentrTPC();

   if (cen < 0 || cen >= 100) { // TPC centrality not defined
      return false;
   }

   mhEvents->Fill(2.5);

   mhVertex->Fill(mPrimaryVertex.Z());
   mhCentrality->Fill(cen);

   return true;
}

float MpdTrackPidMaker::dEdx_sigma_El(float dEdx, float mom) const
{
   if (mom < 0.04) return -999;
   if (dEdx <= 0) return -999;

   dEdx = log(dEdx);

   if (mom < 0.06) mom = 0.06;
   if (mom > 2.0) mom = 2.0;

   float mean[7]  = {6.869878e-003, -1.573146e-001, 1.847371e+000, -9.253678e-001,
                    1.073907e+000, -3.394239e+000, 6.451762e-001};
   float width[7] = {-4.359368e+006, -8.508504e-012, -3.958364e-009, 1.526816e-009,
                     1.353776e-011,  3.426352e-009,  6.591542e-002};

   float mean_exp, width_exp;

   mean_exp =
      mean[0] / mom / mom *
         (mean[1] * log(mom * mom) - mean[2] * mom * mom - mean[3] * mom - mean[4] - mean[5] * mom * mom * mom) +
      mean[6];
   width_exp =
      width[0] / mom / mom *
         (width[1] * log(mom * mom) - width[2] * mom * mom - width[3] * mom - width[4] - width[5] * mom * mom * mom) +
      width[6];

   return (dEdx - mean_exp) / width_exp;
}

float MpdTrackPidMaker::dEdx_sigma_Pi(float dEdx, float mom) const
{
   if (mom < 0.06) return -999;
   if (dEdx <= 0) return -999;

   dEdx = log(dEdx);

   if (mom < 0.07) mom = 0.07;
   if (mom > 3.5) mom = 3.5;

   float mean[7]  = {-3.554770e-002, -2.666590e-001, 5.701508e+000, -5.220967e+000,
                    1.987420e+000,  1.087186e+000,  1.037435e-001};
   float width[7] = {-3.495341e+007, 2.472950e-011,  -1.422238e-010, 3.234257e-010,
                     -1.242149e-010, -5.320861e-011, 7.771059e-002};

   float mean_exp, width_exp;

   mean_exp =
      mean[0] / mom / mom *
         (mean[1] * log(mom * mom) - mean[2] * mom * mom - mean[3] * mom - mean[4] - mean[5] * mom * mom * mom) +
      mean[6];
   width_exp =
      width[0] / mom / mom *
         (width[1] * log(mom * mom) - width[2] * mom * mom - width[3] * mom - width[4] - width[5] * mom * mom * mom) +
      width[6];

   return (dEdx - mean_exp) / width_exp;
}

float MpdTrackPidMaker::dEdx_sigma_K(float dEdx, float mom) const
{
   if (mom < 0.08) return -999;
   if (dEdx <= 0) return -999;

   dEdx = log(dEdx);

   if (mom < 0.09) mom = 0.09;
   if (mom > 3.5) mom = 3.5;

   float mean[7]  = {-4.571263e-004, -1.255796e+001, -7.145174e+002, 1.023127e+003,
                    3.497842e+001,  2.890434e+002,  -3.740436e-003};
   float width[7] = {-9.820512e+006, -1.117400e-010, 1.188118e-009, -1.754282e-009,
                     7.875521e-010,  -4.171753e-010, 8.100840e-002};

   float mean_exp, width_exp;

   mean_exp =
      mean[0] / mom / mom *
         (mean[1] * log(mom * mom) - mean[2] * mom * mom - mean[3] * mom - mean[4] - mean[5] * mom * mom * mom) +
      mean[6];
   width_exp =
      width[0] / mom / mom *
         (width[1] * log(mom * mom) - width[2] * mom * mom - width[3] * mom - width[4] - width[5] * mom * mom * mom) +
      width[6];

   return (dEdx - mean_exp) / width_exp;
}

float MpdTrackPidMaker::dEdx_sigma_P(float dEdx, float mom) const
{
   if (mom < 0.08) return -999;
   if (dEdx <= 0) return -999;

   dEdx = log(dEdx);

   if (mom < 0.09) mom = 0.09;
   if (mom > 3.5) mom = 3.5;

   float mean[7]  = {1.183672e-001, -3.394823e-001, 2.418100e+001, -1.377881e+001,
                    2.534160e+000, -1.563054e+000, 2.010767e+000};
   float width[7] = {-1.488536e+007, -2.075514e-010, 2.418077e-009, -2.996053e-009,
                     1.397806e-009,  -4.161277e-010, 7.174217e-002};

   float mean_exp, width_exp;

   mean_exp =
      mean[0] / mom / mom *
         (mean[1] * log(mom * mom) - mean[2] * mom * mom - mean[3] * mom - mean[4] - mean[5] * mom * mom * mom) +
      mean[6];
   width_exp =
      width[0] / mom / mom *
         (width[1] * log(mom * mom) - width[2] * mom * mom - width[3] * mom - width[4] - width[5] * mom * mom * mom) +
      width[6];

   return (dEdx - mean_exp) / width_exp;
}

float MpdTrackPidMaker::TofdPhiMatch(int charge, float pt, float dphi) const
{
   if (pt < 0.1) return -999;

   float meanphiCoeff[7]  = {2.745900e-001, -1.437619e+000, 2.824655e+000, -2.644983e+000,
                            1.269623e+000, -3.012448e-001, 2.816063e-002};
   float widthphiCoeff[7] = {3.211516e+000, -2.325625e+001, 8.024492e+001, -1.407389e+002,
                             1.320427e+002, -6.304419e+001, 1.202106e+001};

   float meanphi, widthphi;

   if (pt < 3.0) {
      meanphi = meanphiCoeff[0] + meanphiCoeff[1] * pt + meanphiCoeff[2] * pt * pt + meanphiCoeff[3] * pt * pt * pt +
                meanphiCoeff[4] * pt * pt * pt * pt + meanphiCoeff[5] * pt * pt * pt * pt * pt +
                meanphiCoeff[6] * pt * pt * pt * pt * pt * pt;
   } else {
      meanphi = 0.135;
   }

   if (charge < 0) meanphi = -meanphi;

   if (pt < 1.325) {
      widthphi = widthphiCoeff[0] + widthphiCoeff[1] * pt + widthphiCoeff[2] * pt * pt +
                 widthphiCoeff[3] * pt * pt * pt + widthphiCoeff[4] * pt * pt * pt * pt +
                 widthphiCoeff[5] * pt * pt * pt * pt * pt + widthphiCoeff[6] * pt * pt * pt * pt * pt * pt;
   } else {
      widthphi = 0.454 + 0.00815163 * pt;
   }

   return (dphi - meanphi) / widthphi;
}

float MpdTrackPidMaker::TofdZedMatch(float pt, float dz) const
{
   if (pt < 0.1) return -999;

   float widthzCoeff[7] = {1.646438e+000, -9.744003e+000, 3.054315e+001, -4.799731e+001,
                           3.994061e+001, -1.680153e+001, 2.810887e+000};

   float meanz, widthz;

   meanz = 0.0;

   if (pt < 1.5) {
      widthz = widthzCoeff[0] + widthzCoeff[1] * pt + widthzCoeff[2] * pt * pt + widthzCoeff[3] * pt * pt * pt +
               widthzCoeff[4] * pt * pt * pt * pt + widthzCoeff[5] * pt * pt * pt * pt * pt +
               widthzCoeff[6] * pt * pt * pt * pt * pt * pt;
   } else {
      widthz = 0.39844 - 0.00660409 * pt;
   }

   return (dz - meanz) / widthz;
}

float MpdTrackPidMaker::Beta_sigma_El(float beta, float mom) const
{
   if (mom < 0.1) return -999;
   if (mom > 1.5) mom = 1.5;

   float mean[7]  = {3.150000e+003, -5.949618e-009, 5.973763e-002, -4.258060e-007,
                    7.841692e-008, -6.563439e-007, 1.891766e+002};
   float width[7] = {-2.501047e+006, -1.253288e-011, 7.294240e-011, -1.362490e-010,
                     7.314968e-011,  -5.672957e-010, 1.322390e-002};

   float mean_exp, width_exp;

   mean_exp =
      mean[0] / mom / mom *
         (mean[1] * log(mom * mom) - mean[2] * mom * mom - mean[3] * mom - mean[4] - mean[5] * mom * mom * mom) +
      mean[6];
   width_exp =
      width[0] / mom / mom *
         (width[1] * log(mom * mom) - width[2] * mom * mom - width[3] * mom - width[4] - width[5] * mom * mom * mom) +
      width[6];

   return (beta - mean_exp) / width_exp;
}

float MpdTrackPidMaker::Beta_sigma_Pi(float beta, float mom) const
{
   if (mom < 0.12) return -999;
   if (mom > 4.0) mom = 4.0;

   float mean[7]  = {3.150000e+003, -1.248103e-006, 1.108943e+000, -7.829288e-006,
                    7.832307e-006, -4.780232e-007, 3.494167e+003};
   float width[7] = {-3.587374e+007, 1.391129e-012,  3.658606e-011, -2.025494e-011,
                     -7.243473e-013, -3.150429e-011, 1.273604e-002};

   float mean_exp, width_exp;

   mean_exp =
      mean[0] / mom / mom *
         (mean[1] * log(mom * mom) - mean[2] * mom * mom - mean[3] * mom - mean[4] - mean[5] * mom * mom * mom) +
      mean[6];
   width_exp =
      width[0] / mom / mom *
         (width[1] * log(mom * mom) - width[2] * mom * mom - width[3] * mom - width[4] - width[5] * mom * mom * mom) +
      width[6];

   return (beta - mean_exp) / width_exp;
}

float MpdTrackPidMaker::Beta_sigma_K(float beta, float mom) const
{
   if (mom < 0.2) return -999;
   if (mom > 3.25) mom = 3.25;

   float mean[7]  = {3.150000e+003, -3.229103e-006, 1.500830e+000, 5.001743e-005,
                    9.494654e-006, 5.751109e-006,  4.728718e+003};
   float width[7] = {-7.187794e+006, -9.693571e-011, 6.368089e-010, -1.548974e-009,
                     6.073912e-010,  -3.236421e-010, 1.551237e-002};

   float mean_exp, width_exp;

   mean_exp =
      mean[0] / mom / mom *
         (mean[1] * log(mom * mom) - mean[2] * mom * mom - mean[3] * mom - mean[4] - mean[5] * mom * mom * mom) +
      mean[6];
   width_exp =
      width[0] / mom / mom *
         (width[1] * log(mom * mom) - width[2] * mom * mom - width[3] * mom - width[4] - width[5] * mom * mom * mom) +
      width[6];

   return (beta - mean_exp) / width_exp;
}

float MpdTrackPidMaker::Beta_sigma_P(float beta, float mom) const
{
   if (mom < 0.3) return -999;
   if (mom > 4.5) mom = 4.5;

   float mean[7]  = {3.150000e+003,  1.983733e-006, 1.499777e+000, 1.733988e-004,
                    -3.201559e-005, 7.804674e-006, 4.725502e+003};
   float width[7] = {-2.663514e+007, 1.388684e-011,  4.657636e-011, -2.603998e-011,
                     -2.325289e-011, -1.370724e-011, 1.105921e-002};

   float mean_exp, width_exp;

   mean_exp =
      mean[0] / mom / mom *
         (mean[1] * log(mom * mom) - mean[2] * mom * mom - mean[3] * mom - mean[4] - mean[5] * mom * mom * mom) +
      mean[6];
   width_exp =
      width[0] / mom / mom *
         (width[1] * log(mom * mom) - width[2] * mom * mom - width[3] * mom - width[4] - width[5] * mom * mom * mom) +
      width[6];

   return (beta - mean_exp) / width_exp;
}

float MpdTrackPidMaker::EmcdPhiMatch(int charge, float pt, float dphi) const
{
   if (pt < 0.1) return -999;
   if (pt > 1.0) pt = 1.0;

   // dPhi (charge>0)
   float meanphiCoeff[6]  = {-3.231339e-002, 1.078172e-001,  -1.315846e-001,
                            7.660724e-002,  -2.138274e-002, 2.299298e-003};
   float widthphiCoeff[6] = {1.154620e-001,  -7.403730e-001, 2.249572e+000,
                             -3.435599e+000, 2.558980e+000,  -7.380640e-001};

   float meanphi = meanphiCoeff[0] + meanphiCoeff[1] * pt + meanphiCoeff[2] * pt * pt + meanphiCoeff[3] * pt * pt * pt +
                   meanphiCoeff[4] * pt * pt * pt * pt + meanphiCoeff[5] * pt * pt * pt * pt * pt;

   if (charge < 0) meanphi = -meanphi;

   float widthphi = widthphiCoeff[0] + widthphiCoeff[1] * pt + widthphiCoeff[2] * pt * pt +
                    widthphiCoeff[3] * pt * pt * pt + widthphiCoeff[4] * pt * pt * pt * pt +
                    widthphiCoeff[5] * pt * pt * pt * pt * pt;

   return (dphi - meanphi) / widthphi;
}

float MpdTrackPidMaker::EmcdZedMatch(float pt, float dz) const
{
   if (pt < 0.1) return -999;
   if (pt > 1.0) pt = 1.0;

   float widthzCoeff[6] = {8.452819e+000, -2.572417e+001, 4.050086e+001, -3.097804e+001, 1.138158e+001, -1.603740e+000};

   float meanz = 0.;

   float widthz = widthzCoeff[0] + widthzCoeff[1] * pt + widthzCoeff[2] * pt * pt + widthzCoeff[3] * pt * pt * pt +
                  widthzCoeff[4] * pt * pt * pt * pt + widthzCoeff[5] * pt * pt * pt * pt * pt;

   return (dz - meanz) / widthz;
}
