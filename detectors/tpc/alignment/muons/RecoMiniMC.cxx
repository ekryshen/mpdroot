///
/// \brief MiniMC for TPC alignment tests
///
/// \author Valentin Kuzmin, SINP of Moscow State University

#include "TpcMiniMC.h"
#include "RecoMiniMC.h"
#include <iostream>
using namespace std;

#include <TVector2.h>
#include <TVector3.h>
#include <TGeoManager.h>
#include <TGeoPgon.h>
#include <TGeoTube.h>
#include <TMath.h>
#include <TSystem.h>
#include <Riostream.h>
#include <time.h>
#include <TString.h>
#include <TMath.h>
#include <TRotation.h>

void RecoMiniMC::Init(BaseTpcSectorGeo &tpcGeo, Int_t Event1, Int_t Event2, Int_t *Sectors, Int_t MinHits,
                      const char *InDataFile, // full path
                      const char *AInFile,    // full path
                      const char *OutDir,     // full pat
                      const char *Comment, const char *OutDataFile)
{ // only name )
   InHitsFile  = InDataFile;
   OutHitsFile = OutDataFile;
   aInFile     = AInFile;
   outDir      = OutDir;
   event1      = Event1;
   event2      = Event2;
   if (MinHits < 0) {
      newfit  = false;
      minHits = -MinHits;
   } else {
      newfit  = true;
      minHits = MinHits;
   }
   for (Int_t i = 0; i < 24; i++) sectors[i] = Sectors[i];
   comment = Comment;
   if (printL > 0) cout << "RecoMiniMC: input file: " << InHitsFile.Data() << endl;
   fTpcSecGeo = dynamic_cast<TpcSectorGeoAlignmentVAK *>(&tpcGeo);
   if (!fTpcSecGeo) Fatal("RecoMiniMC::Init", " !!! Wrong geometry type !!! ");
   InputMCfile   = new TFile(InHitsFile.Data(), "READ");
   tin_alignment = (TTree *)InputMCfile->Get("alignment");
   tin_alignment->SetBranchAddress("R0_A", &old_R0_A);
   tin_alignment->SetBranchAddress("alpha_A", &old_alpha_A);
   tin_alignment->SetBranchAddress("beta_A", &old_beta_A);
   tin_alignment->SetBranchAddress("gamma_A", &old_gamma_A);
   tin_alignment->GetEvent(0);
   delete tin_alignment;

   // initialize the TPC Sectot geometry with the input file alignment
   fTpcSecGeo->TpcSectorGeoA(old_R0_A, old_alpha_A, old_beta_A, old_gamma_A);
   for (Int_t isec = 0; isec < 24; isec++) { // save Alignment used for the sample
      fTpcSecGeo->SectorBackTransformation(isec, oldG0_gsc[isec], old_glo2loc[isec]);
      fTpcSecGeo->SectorTransformation(isec, oldR0_shift[isec], old_loc2glo[isec]);
   }

   // read new alignment to recalculate global hita
   if (printL > 0) cout << "RecoMiniMC: TPC alignment from: " << aInFile.Data() << endl;
   alignmentF    = new TFile(aInFile, "READ");
   tin_alignment = (TTree *)alignmentF->Get("alignment");
   tin_alignment->SetBranchAddress("R0_A", &R0_A);
   tin_alignment->SetBranchAddress("alpha_A", &alpha_A);
   tin_alignment->SetBranchAddress("beta_A", &beta_A);
   tin_alignment->SetBranchAddress("gamma_A", &gamma_A);
   tin_alignment->GetEvent(0);
   delete tin_alignment;
   fTpcSecGeo->TpcSectorGeoA(R0_A, alpha_A, beta_A, gamma_A);
   for (Int_t isec = 0; isec < 24; isec++) {
      fTpcSecGeo->SectorBackTransformation(isec, newG0_gsc[isec], new_glo2loc[isec]);
      fTpcSecGeo->SectorTransformation(isec, newR0_shift[isec], new_loc2glo[isec]);
   }
   alignmentF->Close();
   delete alignmentF;

   if (printL > 0) {
      printf("RecoMiniMC:================== input data ======================\n");
      // fullChi2 TRee
      TTree *tin_fullChi2;
      tin_fullChi2 = (TTree *)InputMCfile->Get("Chi2");
      tin_fullChi2->SetBranchAddress("fNhits", &fNhits);
      tin_fullChi2->SetBranchAddress("fSumChi2", &fSumChi2);
      tin_fullChi2->SetBranchAddress("secNhits", &fNhitSec);
      tin_fullChi2->SetBranchAddress("secSumChi2", &fSecChi2);
      tin_fullChi2->GetEvent(0);
      Double_t sum = 0;
      Long64_t nhs = 0;
      for (Int_t s = 0; s < 24; s++)
         if (Sectors[s]) {
            sum += fSecChi2[s];
            nhs += fNhitSec[s];
            fSecChi2[s] = (fNhitSec[s] > 0) ? fSecChi2[s] / fNhitSec[s] : 0;
            printf("sector=%2d\n", s);
            printf("old R0(%7.4f,%7.4f,%7.4f)  abg(%7.4f,%7.4f,%7.4f)\n", old_R0_A[0][s], old_R0_A[1][s],
                   old_R0_A[2][s], old_alpha_A[s], old_beta_A[s], old_gamma_A[s]);
            printf("old chi2=%18.15f   nhits=%lld\n", fSecChi2[s], fNhitSec[s]);
            printf("new R0(%7.4f,%7.4f,%7.4f)  abg(%7.4f,%7.4f,%7.4f)\n", R0_A[0][s], R0_A[1][s], R0_A[2][s],
                   alpha_A[s], beta_A[s], gamma_A[s]);
         }
      sum = (sum > 0) ? sum / nhs : 0;
      printf("======= old chi2/ndf sum %18.15f  nhits=%lld =========\n", sum, nhs);
   }

   for (Int_t is = 0; is < 24; is++) {
      fNhitSec[is] = 0;
      fSecChi2[is] = 0;
   }
   fNhits    = 0;
   fSumChi2  = 0;
   proEvents = 0;
}

//______________________________________________________________________
// simulate fnEvents and write data to OutPut file
void RecoMiniMC::reco()
{
   // create the output file name
   stringstream ss;
   if (!outDir.EqualTo(" ")) { // only incomplate output file name
      Ssiz_t  ld    = InHitsFile.Last('/');
      Ssiz_t  ln    = InHitsFile.Last('.');
      TString fName = InHitsFile(ld + 1, ln - ld - 1);
      if (!OutHitsFile.EqualTo(" ")) {
         ss << outDir.Data() << "/" << OutHitsFile.Data();
      } else {
         ss << outDir.Data() << "/Reco_" << fName.Data() << comment;
         if (!aInFile.EqualTo(" ")) { // a special alignment
            ld    = aInFile.Last('/');
            ln    = aInFile.Last('.');
            fName = aInFile(ld + 1, ln - ld - 1);
            ss << fName.Data();
         }
      }
      ss << ".root";
   } else {
      ss << OutHitsFile.Data();
   }
   // write the alignment to the output file
   TString *outFile = new TString(ss.str());
   if (printL > 0) {
      cout << "RecoMiniMC::reco:   inputfile: " << InHitsFile.Data() << endl;
      cout << "RecoMiniMC::reco:  outputfile: " << outFile->Data() << endl;
      cout << "RecoMiniMC::reco: InAlignment: " << aInFile.Data() << endl;
   }
   OutputReco = new TFile(outFile->Data(), "RECREATE");
   delete outFile;
   // wtite old alignment to the output file
   old_alignment = new TTree("old_alignment", "alignment_parameters");
   old_alignment->Branch("R0_A", &old_R0_A, "old_R0_A[3][24]/D");
   old_alignment->Branch("alpha_A", &old_alpha_A, "old_alpha_A[24]/D");
   old_alignment->Branch("beta_A", &old_beta_A, "old_beta_A[24]/D");
   old_alignment->Branch("gamma_A", &old_gamma_A, "old_gamma_A[24]/D");
   old_alignment->Fill();
   old_alignment->Write();
   delete old_alignment;

   // wtite new alignment to the output file
   tout_alignment = new TTree("alignment", "alignment_parameters");
   tout_alignment->Branch("R0_A", &R0_A, "R0_A[3][24]/D");
   tout_alignment->Branch("alpha_A", &alpha_A, "alpha_A[24]/D");
   tout_alignment->Branch("beta_A", &beta_A, "beta_A[24]/D");
   tout_alignment->Branch("gamma_A", &gamma_A, "gamma_A[24]/D");
   tout_alignment->Fill();
   tout_alignment->Write();
   delete tout_alignment;
   // for(Int_t i=0;i<12;i++) printf("RECO: sec=%d outputA(%f %f %f  %f %f %f)\n",
   // i,R0_A[0][i],R0_A[1][i],R0_A[2][i],alpha_A[i],beta_A[i],gamma_A[i]);
   //  read parametrs from the input file
   InputMCfile->cd();
   tin_parameters = (TTree *)InputMCfile->Get("parameters");

   tin_parameters->SetBranchAddress("b_mf", &b_mf);
   tin_parameters->SetBranchAddress("e_ef", &e_ef);
   tin_parameters->SetBranchAddress("r0_beam", &r0_beam);
   tin_parameters->SetBranchAddress("p_las", &p_las);
   tin_parameters->SetBranchAddress("d_track", &d_track);
   tin_parameters->SetBranchAddress("p_min", &p_min);
   tin_parameters->SetBranchAddress("p_max", &p_max);
   tin_parameters->SetBranchAddress("min_padS", &min_padS);
   tin_parameters->SetBranchAddress("max_padS", &max_padS);
   tin_parameters->SetBranchAddress("sig_padS", &sig_padS);
   tin_parameters->SetBranchAddress("sig_timeBins", &sig_timeBins);
   tin_parameters->SetBranchAddress("sig_Z_pos", &sig_Z_pos);
   tin_parameters->SetBranchAddress("sector0", &sector0);
   tin_parameters->SetBranchAddress("sector1", &sector1);
   tin_parameters->SetBranchAddress("k_par", &k_par);
   tin_parameters->SetBranchAddress("maxPcluster", &maxPcluster);
   tin_parameters->SetBranchAddress("fSeed", &fSeed);
   tin_parameters->SetBranchAddress("beam_point", &beam_point);
   tin_parameters->SetBranchAddress("MFieldMode", &MFieldMode);
   tin_parameters->GetEvent(0);
   delete tin_parameters;

   if (MFieldMode) {
      bMf = TVector3(b_mf[0], b_mf[1], b_mf[2]);
      Double_t phi_B;
      if (bMf.X() * bMf.X() + bMf.Y() * bMf.Y() == 0)
         phi_B = 0;
      else
         phi_B = TMath::ATan2(bMf.Y(), bMf.X());
      B2G.SetToIdentity();
      B2G.RotateZ(phi_B);
      Double_t theta_B = TMath::ACos(bMf.Z() / bMf.Mag());
      B2G.RotateY(theta_B);
      G2B = B2G.Inverse();
   }

   eEf    = TVector3(e_ef[0], e_ef[1], e_ef[2]);
   r0Beam = TVector3(r0_beam[0], r0_beam[1], r0_beam[2]);
   pLas   = TVector3(p_las[0], p_las[1], p_las[2]);

   // write parametrs to the output file
   OutputReco->cd();
   tout_parameters = new TTree("parameters", "miniMC_parameters");

   tout_parameters->Branch("b_mf", &b_mf, "b_mf[3]/D");
   tout_parameters->Branch("e_ef", &e_ef, "e_ef[3]/D");
   tout_parameters->Branch("r0_beam", &r0_beam, "r0_beam[3]/D");
   tout_parameters->Branch("p_las", &p_las, "p_las[3]/D");
   tout_parameters->Branch("d_track", &d_track, "d_track/D");
   tout_parameters->Branch("p_min", &p_min, "p_min/D");
   tout_parameters->Branch("p_max", &p_max, "p_max/D");
   tout_parameters->Branch("min_padS", &min_padS, "min_padS[2]/D");
   tout_parameters->Branch("max_padS", &max_padS, "max_padS[2]/D");
   tout_parameters->Branch("sig_padS", &sig_padS, "sig_padS[2]/D");
   tout_parameters->Branch("sig_timeBins", &sig_timeBins, "sig_timeBins/D");
   tout_parameters->Branch("sig_Z_pos", &sig_Z_pos, "sig_Z_pos/D");
   tout_parameters->Branch("sector0", &sector0, "sector0/I");
   tout_parameters->Branch("sector1", &sector1, "sector1/I");
   tout_parameters->Branch("k_par", &k_par, "k_par/I");
   tout_parameters->Branch("maxPcluster", &maxPcluster, "maxPcluster/I");
   tout_parameters->Branch("fSeed", &fSeed, "fSeed/i");
   tout_parameters->Branch("beam_point", &beam_point, "beam_point/O");
   tout_parameters->Branch("MFieldMode", &MFieldMode, "MFieldMode/O");
   // write parametrs to the output file
   tout_parameters->Fill();
   tout_parameters->Write();
   delete tout_parameters;

   // data TRee information from the input file
   tin_data = (TTree *)InputMCfile->Get("data");

   tin_data->SetBranchAddress("trackPin", &trackPin);
   tin_data->SetBranchAddress("trackPout", &trackPout);
   tin_data->SetBranchAddress("Chi2Fit", &Chi2Fit);
   tin_data->SetBranchAddress("Par0R", &trackDiff[0]);
   tin_data->SetBranchAddress("Par1R", &trackDiff[1]);
   tin_data->SetBranchAddress("Par2R", &trackDiff[2]);
   tin_data->SetBranchAddress("Par3R", &trackDiff[3]);
   tin_data->SetBranchAddress("Par4R", &trackDiff[4]);
   tin_data->SetBranchAddress("Par5R", &trackDiff[5]);

   tin_data->SetBranchAddress("MpdTpcMiniMC_Nhits", &MpdTpcMiniMC_Nhits);
   tin_data->SetBranchAddress("MpdTpcMiniMC_iSec", &MpdTpcMiniMC_iSec);
   tin_data->SetBranchAddress("MpdTpcMiniMC_XHits", &MpdTpcMiniMC_XHits);
   tin_data->SetBranchAddress("MpdTpcMiniMC_YHits", &MpdTpcMiniMC_YHits);
   tin_data->SetBranchAddress("MpdTpcMiniMC_ZHits", &MpdTpcMiniMC_ZHits);
   tin_data->SetBranchAddress("MpdTpcMiniMC_WHits", &MpdTpcMiniMC_WHits);
   tin_data->SetBranchAddress("MpdTpcMiniMC_chi2s", &MpdTpcMiniMC_chi2s);

   // data TRee information from the input file
   tout_data = new TTree("data", "recominiMC_data");

   tout_data->Branch("trackPin", &trackPin, "trackPin[6]/D");
   tout_data->Branch("trackPout", &trackPout, "trackPout[6]/D");
   tout_data->Branch("Chi2Fit", &Chi2Fit, "Chi2Fit/D");
   tout_data->Branch("Par0R", &trackDiff[0], "trackDiff/D");
   tout_data->Branch("Par1R", &trackDiff[1], "trackDiff/D");
   tout_data->Branch("Par2R", &trackDiff[2], "trackDiff/D");
   tout_data->Branch("Par3R", &trackDiff[3], "trackDiff/D");
   tout_data->Branch("Par4R", &trackDiff[4], "trackDiff/D");
   tout_data->Branch("Par5R", &trackDiff[5], "trackDiff/D");

   tout_data->Branch("MpdTpcMiniMC_Nhits", &RecoMiniMC_Nhits, "RecoMiniMC_Nhits/I");
   tout_data->Branch("MpdTpcMiniMC_iSec", &RecoMiniMC_iSec, "RecoMiniMC_iSec[RecoMiniMC_Nhits]/I");
   tout_data->Branch("MpdTpcMiniMC_XHits", &RecoMiniMC_XHits, "RecoMiniMC_XHits[RecoMiniMC_Nhits]/D");
   tout_data->Branch("MpdTpcMiniMC_YHits", &RecoMiniMC_YHits, "RecoMiniMC_YHits[RecoMiniMC_Nhits]/D");
   tout_data->Branch("MpdTpcMiniMC_ZHits", &RecoMiniMC_ZHits, "RecoMiniMC_ZHits[RecoMiniMC_Nhits]/D");
   tout_data->Branch("MpdTpcMiniMC_WHits", &RecoMiniMC_WHits, "RecoMiniMC_WHits[RecoMiniMC_Nhits]/D");
   tout_data->Branch("MpdTpcMiniMC_chi2s", &RecoMiniMC_chi2s, "RecoMiniMC_chi2s[RecoMiniMC_Nhits]/D");

   // fullChi2 TRee
   tout_fullChi2 = new TTree("Chi2", "fullChi2");

   tout_fullChi2->Branch("fNhits", &fNhits, "fNhits/L");
   tout_fullChi2->Branch("fSumChi2", &fSumChi2, "fSumChi2/D");
   tout_fullChi2->Branch("secNhits", &fNhitSec, "fNhitSec[24]/L");
   tout_fullChi2->Branch("secSumChi2", &fSecChi2, "fSecChi2[24]/D");

   // ****************  cicle over events ********************
   if (event2 == 0) { // if not the constructor #2
      event1 = 0;
      event2 = tin_data->GetEntries();
   }
   if (printL > 0) cout << "RecoMiniMC::reco: analyze events: " << event1 << "-" << event2 << endl;

   InputMCfile->cd();
   for (Int_t iEvt = 0; iEvt < event1; iEvt++) tin_data->GetEvent(iEvt);
   Nreconstructed = 0;
   for (Int_t iEvt = event1; iEvt < event2; iEvt++) {
      RecoMiniMC_IEVT = iEvt;
      if (printL > 2)
         if (iEvt % printL == 0 && iEvt != 0) {
            cout << "RecoMiniMC: event=" << iEvt << ", processed=" << proEvents + 1 << endl;
         }
      InputMCfile->cd();
      tin_data->GetEvent(iEvt);
      if (MpdTpcMiniMC_Nhits >= minHits) {
         proEvents++;
         convert();
         if (!newfit) {
            //        for(Int_t i=0;i<6;i++) trackPout[i]=trackPin[i];
            if (!MFieldMode)
               putChi2L();
            else
               putChi2H();
         } else {
            if (!MFieldMode) {
               if (fitLine() != 0) continue;
               putChi2L();
            } else {
               if (fitHelix() != 0) continue;
               putChi2H();
            }
         }
         Nreconstructed++;
         OutputReco->cd();
         tout_data->Fill();
      }
   }
   OutputReco->cd();
   tout_data->Write();
   if (printL > 0)
      printf("*******  %d events reconstucted of %d produced (%%%4.1f) **********\n", Nreconstructed, proEvents,
             float_t(Nreconstructed) / float_t(proEvents) * 100);

   fNhits   = 0;
   fSumChi2 = 0;
   for (Int_t is = 0; is < 24; is++) {
      if (sectors[is]) {
         fNhits += fNhitSec[is];
         fSumChi2 += fSecChi2[is];
      }
   }

   if (fNhits > 0) fSumChi2 /= fNhits;
   tout_fullChi2->Fill();
   tout_fullChi2->Write();

   //    OutputReco->Close();
   delete tout_data;
   delete OutputReco;
   //    InputMCfile->Close();
   delete tin_data;
   delete InputMCfile;
   Long64_t fHts = fNhits;
   if (printL > 0)
      cout << "RecoMiniMC: ========== Full Chi2 = " << fSumChi2 << "  for " << fHts << " hits in " << proEvents
           << " events =========" << endl;
}

//______________________________________________________________________
// Chi2 for the clusters fit by the direct line
double RecoMiniMC_Chi2Line(const double *par)
{
   double chisq = 0;
   double delta;
   Int_t  nhits = RecoMiniMC_Nhits;
   for (int i = 0; i < nhits; i++) {
      // cout<<"chi2line RecoMiniMC_iSec["<<i<<"]="<<RecoMiniMC_iSec[i]<<endl;
      if (sectors[RecoMiniMC_iSec[i]]) {
         double X = RecoMiniMC_XHits[i];
         double Y = RecoMiniMC_YHits[i];
         double Z = RecoMiniMC_ZHits[i];
         double W = RecoMiniMC_WHits[i];

         double ex = TMath::Sin(par[3]) * TMath::Cos(par[4]);
         double ey = TMath::Sin(par[3]) * TMath::Sin(par[4]);
         double ez = TMath::Cos(par[3]);
         double t  = ex * (X - par[0]) + ey * (Y - par[1]) + ez * (Z - par[2]);

         double X1 = par[0] + ex * t;
         double Y1 = par[1] + ey * t;
         double Z1 = par[2] + ez * t;

         double dx = X - X1;
         double dy = Y - Y1;
         double dz = Z - Z1;
         delta     = TMath::Sqrt(dx * dx + dy * dy + dz * dz) / W;
         chisq += delta * delta;
         // cout<<i<<"chi2line (X,Y,Z) ("<<X<<" "<<Y<<" "<<Z<<")"<<" W="<<W<<endl;
      }
   }
   // cout<<par[0]<<" "<<par[1]<<" "<<par[2]<<" "<<par[3]<<" "<<par[4]<<" Chi2Line="<<chisq<<endl;
   // cout<<chisq<<endl;
   return chisq;
}

//______________________________________________________________________
// find the fit by the direct line
Int_t RecoMiniMC::fitLine()
{
   // Choose method upon creation between:
   // kMigrad, kSimplex, kCombined,
   // kScan, kFumili
   static ROOT::Minuit2::Minuit2Minimizer minFit(ROOT::Minuit2::kMigrad);

   ROOT::Math::Functor f(&RecoMiniMC_Chi2Line, 5);
   double              step[5] = {0.1, 0.1, 0.1, 0.01, 0.01};
   double              limL[5] = {-133, -133, -170, 0, -pi};
   double              limU[5] = {133, 133, 170, pi, pi};
   double              par[5];
   Int_t               imin, imax, iminx, imaxx, iminy, imaxy, iminz, imaxz;
   Double_t            xmn = 200, xmx = -200, ymn = 200, ymx = -200, zmn = 200, zmx = -200, X = 0, Y = 0, Z = 0;
   // printf("fitLine: nHits=%d \n",RecoMiniMC_Nhits);
   for (Int_t i = 0; i < RecoMiniMC_Nhits; i++) {
      // printf("fitLine: i=%d   r(%f,%f,%f) sector=%d\n",i,
      // RecoMiniMC_XHits[i],RecoMiniMC_YHits[i],RecoMiniMC_ZHits[i],MpdTpcMiniMC_iSec[i]);
      if (RecoMiniMC_ZHits[i] > zmx) {
         zmx   = RecoMiniMC_ZHits[i];
         imaxz = i;
      }
      if (RecoMiniMC_ZHits[i] < zmn) {
         zmn   = RecoMiniMC_ZHits[i];
         iminz = i;
      }

      if (RecoMiniMC_XHits[i] > xmx) {
         xmx   = RecoMiniMC_XHits[i];
         imaxx = i;
      }
      if (RecoMiniMC_XHits[i] < xmn) {
         xmn   = RecoMiniMC_XHits[i];
         iminx = i;
      }

      if (RecoMiniMC_YHits[i] > ymx) {
         ymx   = RecoMiniMC_YHits[i];
         imaxy = i;
      }
      if (RecoMiniMC_YHits[i] < ymn) {
         ymn   = RecoMiniMC_YHits[i];
         iminy = i;
      }

      X += RecoMiniMC_XHits[i];
      Y += RecoMiniMC_YHits[i];
      Z += RecoMiniMC_ZHits[i];
   }

   Double_t dx = TMath::Abs(RecoMiniMC_XHits[imaxx] - RecoMiniMC_XHits[iminx]);
   Double_t dy = TMath::Abs(RecoMiniMC_YHits[imaxy] - RecoMiniMC_YHits[iminy]);
   Double_t dz = TMath::Abs(RecoMiniMC_ZHits[imaxz] - RecoMiniMC_ZHits[iminz]);
   if (dz > dx && dz > dy) {
      imin = iminz;
      imax = imaxz;
   } else {
      if (dy > dz && dy > dx) {
         imin = iminy;
         imax = imaxy;
      } else {
         imin = iminx;
         imax = imaxx;
      }
   }

   X /= RecoMiniMC_Nhits;
   Y /= RecoMiniMC_Nhits;
   Z /= RecoMiniMC_Nhits;
   TVector3 e =
      TVector3(RecoMiniMC_XHits[imax] - RecoMiniMC_XHits[imin], RecoMiniMC_YHits[imax] - RecoMiniMC_YHits[imin],
               RecoMiniMC_ZHits[imax] - RecoMiniMC_ZHits[imin]);

   if (e.Y() > 0) e = -e;
   e           = e.Unit();
   trackPin[0] = X;
   trackPin[1] = Y;
   trackPin[2] = Z;
   trackPin[3] = TMath::ACos(e.Z());
   trackPin[4] = TMath::ATan2(e.Y(), e.X());
   par[0]      = trackPin[0];
   par[1]      = trackPin[1];
   par[2]      = trackPin[2];
   par[3]      = trackPin[3];
   par[4]      = trackPin[4];
   trackPin[5] = 0;
   for (Int_t i = 0; i < 5; i++) {
      if (i < 3) {
         limL[i] = par[i] - 15;
         limU[i] = par[i] + 15;
      } else {
         limL[i] = par[i] - 0.05;
         limU[i] = par[i] + 0.05;
      }
   };

   minFit.SetFunction(f);

   // Set the free variables to be minimized!
   minFit.SetLimitedVariable(0, "x0", par[0], step[0], limL[0], limU[0]);
   minFit.SetLimitedVariable(1, "y0", par[1], step[1], limL[1], limU[1]);
   minFit.SetLimitedVariable(2, "z0", par[2], step[2], limL[2], limU[2]);
   minFit.SetLimitedVariable(3, "theta", par[4], step[3], limL[3], limU[3]);
   minFit.SetLimitedVariable(4, "phi", par[4], step[4], limL[4], limU[4]);

   Bool_t res = minFit.Minimize();
   //      if(minFit.MinValue()/RecoMiniMC_Nhits>EvalpHit) {
   //        minFit.PrintResults();
   //        return 1;
   //      }
   //    minFit.PrintResults();
   const double *xs = minFit.X();
   trackPin[0]      = trackPout[0];
   trackPout[1]     = trackPin[0];
   trackPout[2]     = trackPin[0];
   trackPout[3]     = trackPin[0];
   trackPout[4]     = trackPin[0];

   trackPout[0] = xs[0];
   trackPout[1] = xs[1];
   trackPout[2] = xs[2];
   trackPout[3] = xs[3];
   trackPout[4] = xs[4];
   trackPout[5] = minFit.NCalls();
   // printf("event=%d Par(%9.4f %9.4f %9.4f %9.4f %9.4f)\n",Nreconstructed,
   // trackPout[0],trackPout[1],trackPout[2],r2d*trackPout[3],r2d*trackPout[4]);
   for (Int_t k = 0; k < 6; k++) trackDiff[k] = trackPout[k] - trackPin[k];
   return 0;
   //    }
   //    else return 1;
}

//______________________________________________________________________
// put chi2 values
void RecoMiniMC::putChi2L()
{
   Double_t chisq = 0;
   Double_t delta;
   TVector3 x1;
   Int_t    nhits = RecoMiniMC_Nhits;
   Double_t ex    = TMath::Sin(trackPout[3]) * TMath::Cos(trackPout[4]);
   Double_t ey    = TMath::Sin(trackPout[3]) * TMath::Sin(trackPout[4]);
   Double_t ez    = TMath::Cos(trackPout[3]);
   // printf("tPout(%f,%f,%f,%f,%f)\n",trackPout[0],trackPout[1],trackPout[2],trackPout[3],trackPout[4]);
   TVector3 e  = TVector3(ex, ey, ez);
   TVector3 x0 = TVector3(trackPout[0], trackPout[1], trackPout[2]);
   for (int i = 0; i < nhits; i++) {
      Double_t X          = RecoMiniMC_XHits[i];
      Double_t Y          = RecoMiniMC_YHits[i];
      Double_t Z          = RecoMiniMC_ZHits[i];
      Double_t W          = RecoMiniMC_WHits[i];
      x1                  = TVector3(RecoMiniMC_XHits[i], RecoMiniMC_YHits[i], RecoMiniMC_ZHits[i]);
      x1                  = e.Cross(x0 - x1);
      delta               = x1.Mag() / W;
      RecoMiniMC_chi2s[i] = delta * delta;
      chisq += delta * delta;
      // printf("putChi2: RecoMiniMC_iSec[%d] delta^2=%f  chi2=%f\n",i,delta*delta,chisq);
      fNhitSec[RecoMiniMC_iSec[i]]++;
      fSecChi2[RecoMiniMC_iSec[i]] += delta * delta;
   }
   // printf("putChi2: oldChi2Fit=%f newChi2Fit=%f dif=%7.4f\n",
   // Chi2Fit,chisq/RecoMiniMC_Nhits,Chi2Fit-chisq/RecoMiniMC_Nhits);
   Chi2Fit = chisq / RecoMiniMC_Nhits;
   // printf("putChi2: Chi2Fit=%f\n",Chi2Fit);
}

//______________________________________________________________________
// Chi2 for the circle projection
Double_t RecoMiniMC_Chi2Circle(const Double_t *par)
{
   //    par[0] - Xc; par[1] - Yc; par[2] - R; par[3] - phi0
   double t;
   double xc    = par[0];
   double yc    = par[1];
   double R     = par[2];
   double phi0  = par[3];
   double X     = xc + R * TMath::Cos(phi0);
   double Y     = yc + R * TMath::Sin(phi0);
   double delta = TMath::Sqrt(X * X + Y * Y) / 0.2;
   double chisq = delta * delta;
   for (int i = 0; i < RecoMiniMC_Nhits; i++) {
      X          = RecoMiniMC_XBHits[i];
      Y          = RecoMiniMC_YBHits[i];
      double dfi = phi0 - TMath::ATan2(Y - yc, X - xc);
      if (RecoMiniMC_qHelix > 0) {
         if (dfi < 0)
            t = TMath::TwoPi() + dfi;
         else
            t = dfi;
      } else {
         if (dfi > 0)
            t = TMath::TwoPi() - dfi;
         else
            t = -dfi;
      }
      double dx = X - (xc + R * TMath::Cos(phi0 - RecoMiniMC_qHelix * t));
      double dy = Y - (yc + R * TMath::Sin(phi0 - RecoMiniMC_qHelix * t));
      // printf("i=%2d t=%6.1f tF=%6.1f H(%6.1f %6.1f) T(%6.1f %6.1f) dx=%f dy=%f chi2=%f\n",i,
      // t*57.29577951,(phi0-RecoMiniMC_qHelix*t)*57.29577951,X,Y,X-dx,Y-dy,dx,dy,chisq);
      delta = TMath::Sqrt(dx * dx + dy * dy) / RecoMiniMC_WHits[i];
      chisq += delta * delta;
   }
   // printf("xc,yc(%6.1f %6.1f) R=%6.1f phi0=%6.1f chi2=%f\n",
   // xc,yc,R,phi0*57.29577951,chisq);
   return chisq;
}

//______________________________________________________________________
// Chi2 for the clusters fit by the helix
Double_t RecoMiniMC_zHelix(const Double_t *par)
{

   double chisq = 0;
   double xc    = trackParOut[0];
   double yc    = trackParOut[1];
   double R     = TMath::Abs(trackParOut[3]);
   double phi0  = trackParOut[4];
   // printf("xc,yc(%6.1f %6.1f) R=%6.1f phi0=%6.1f\n",xc,yc,R,phi0*57.29577951);
   for (int i = 0; i < RecoMiniMC_Nhits; i++) {
      double t     = RecoMiniMC_tHelix[i];
      double dx    = RecoMiniMC_XBHits[i] - (xc + R * TMath::Cos(phi0 - RecoMiniMC_qHelix * t));
      double dy    = RecoMiniMC_YBHits[i] - (yc + R * TMath::Sin(phi0 - RecoMiniMC_qHelix * t));
      double dz    = (RecoMiniMC_ZBHits[i] - (par[0] + par[1] * RecoMiniMC_tHelix[i]));
      double delta = TMath::Sqrt(dx * dx + dy * dy + dz * dz) / RecoMiniMC_WHits[i];
      // printf("i=%2d t=%6.1f (%8.1f %8.1f) dx=%f dy=%f dz=%f chi2=%f\n",i,
      // RecoMiniMC_tHelix[i],par[0],par[1],dx,dy,dz,chisq);
      chisq += delta * delta;
   }
   // printf("par(%lf %lf) chi2=%lf\n",par[0],par[1],chisq);
   return chisq;
}

//______________________________________________________________________
// Reconstruction of the helix parameters
Int_t RecoMiniMC::fitHelix()
{
   /*
     in the BCS the helix equation:
       x=p[0]+p[3]*cos(p[4]+q*t)
       y=p[1]+p[3]*sin(p[4]+q*t)
       z=p[2]+p[5]*t
   */
   Int_t    imin, imax;
   Double_t zmin = 200, zmax = 0;
   for (Int_t i = 0; i < RecoMiniMC_Nhits; i++) {
      if (RecoMiniMC_ZBHits[i] < zmin) {
         imin = i;
         zmin = RecoMiniMC_ZBHits[i];
      }
      if (RecoMiniMC_ZBHits[i] > zmax) {
         imax = i;
         zmax = RecoMiniMC_ZBHits[i];
      }
      // printf("i=%2d sec=%2d XYZ(%8.2f %8.2f %8.2f)\n",i,RecoMiniMC_iSec[i],
      // RecoMiniMC_XBHits[i],RecoMiniMC_YBHits[i],RecoMiniMC_ZBHits[i]);
   }
   if (zmin < 0) {
      Int_t iw    = imin;
      imin        = imax;
      imax        = iw;
      Double_t zw = zmin;
      zmin        = zmax;
      zmax        = zw;
   }
   Double_t x1, y1, z1, x2, y2, z2, c2, x3, y3, z3, xw, yw, t, t1, t2, xc = 0, yc = 0;
   x1 = 0;
   y1 = 0;
   x2 = RecoMiniMC_XBHits[imin];
   y2 = RecoMiniMC_YBHits[imin];
   z2 = RecoMiniMC_ZBHits[imin];
   x3 = RecoMiniMC_XBHits[imax];
   y3 = RecoMiniMC_YBHits[imax];
   z3 = RecoMiniMC_ZBHits[imax];

   CircleCenter(x1, y1, x2, y2, x3, y3, xc, yc),
      // printf("x1,y1,x2,y2,x3,y3,xc,yc(%7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f\n",
      // x1,y1,x2,y2,x3,y3,xc,yc);
      //     xc+=xw; yc+=yw;
      //     CircleCenter(x2,y2,x3,y3,x1,y1,xw,yw),
      //     xc+=xw; yc+=yw;
      //     CircleCenter(x3,y3,x1,y1,x2,y2,xw,yw),
      //     xc+=xw; yc+=yw;
      //     xc/=3; yc/=3;*/

      // determ the particle sign
      t1 = TMath::ATan2(y2 - yc, x2 - xc);
   t2    = TMath::ATan2(y3 - yc, x3 - xc);
   if (t1 > t2) {
      if (t1 - t2 > pi) {
         RecoMiniMC_qHelix = -1;
         t                 = twopi - (t1 - t2);
      } else {
         RecoMiniMC_qHelix = 1;
         t                 = t1 - t2;
      }
   } else {
      if (t2 - t1 > pi) {
         RecoMiniMC_qHelix = 1;
         t                 = twopi - (t2 - t1);
      } else {
         RecoMiniMC_qHelix = -1;
         t                 = t2 - t1;
      }
   }
   // printf("fitHelix: R1(%8.1f,%8.1f) R2(%8.1f %8.1f)  q=%d\n",x2,y2,x3,y3,RecoMiniMC_qHelix);

   // Choose method upon creation between:
   // kMigrad, kSimplex, kCombined,
   // kScan, kFumili
   static ROOT::Minuit2::Minuit2Minimizer minFitc(ROOT::Minuit2::kMigrad);
   //  ================ find 4 circle parameters =================
   ROOT::Math::Functor fcn1(&RecoMiniMC_Chi2Circle, 4);
   //                     Xc    Yc    R     Phi0
   double step[6] = {1, 1, 1, 0.02, 0, 0};
   double limU[6] = {2000, 2000, 2000, TMath::Pi(), 0, 0};
   double limL[6] = {-2000, -2000, -2000, -TMath::Pi(), 0, 0};
   double par[6];

   par[0] = xc;                             //
   par[1] = yc;                             //  - coordinates of the muon start point
   par[2] = TMath::Sqrt(xc * xc + yc * yc); // radius of the circle in XY plane of BCS
   par[3] = TMath::ATan2(-yc, -xc);         // initial muon angle (XY plane of BCS)

   // printf("PinHelix: par(%8.1f %8.1f %8.1f %6.1f)  q=%d nhits=%d\n",
   // trackPin[0],trackPin[1],trackPin[3],r2d*trackPin[4],
   // RecoMiniMC_qHelix,RecoMiniMC_Nhits);
   // printf("fitHelix: par(%8.1f %8.1f %8.1f %6.1f)  q=%d\n",
   // par[0],par[1],par[2],r2d*par[3],RecoMiniMC_qHelix);

   minFitc.SetFunction(fcn1);

   // Set the free variables to be minimized!
   minFitc.SetLimitedVariable(0, "x0", par[0], step[0], limL[0], limU[0]);
   minFitc.SetLimitedVariable(1, "y0", par[1], step[1], limL[1], limU[1]);
   minFitc.SetLimitedVariable(2, "R", par[2], step[2], limL[2], limU[2]);
   minFitc.SetLimitedVariable(3, "phi", par[3], step[3], limL[3], limU[3]);

   Bool_t        res = minFitc.Minimize();
   const double *xs1 = minFitc.X();
   // printf("outHelix: par(%8.1f %8.1f %8.1f %6.1f)  q=%d\n",
   // xs1[0],xs1[1],xs1[2],r2d*xs1[3],RecoMiniMC_qHelix);
   // printf("chi2min4=%f\n",minFitc.MinValue()/RecoMiniMC_Nhits);
   // cout<<"res="<<res<<endl<<endl;
   if (minFitc.MinValue() / RecoMiniMC_Nhits > EvalpHit) return 1;
   //      minFitc.PrintResults();
   xc = trackParOut[0] = trackPout[0] = xs1[0];
   yc = trackParOut[1] = trackPout[1] = xs1[1];
   trackPout[3] = trackParOut[3] = RecoMiniMC_qHelix * xs1[2];
   Double_t phi0 = trackParOut[4] = trackPout[4] = xs1[3];

   trackDiff[0] = trackPout[0] - trackPin[0];
   trackDiff[1] = trackPout[1] - trackPin[1];
   trackDiff[3] = trackPout[3] - trackPin[3];
   trackDiff[4] = trackPout[4] - trackPin[4];

   //  ================ find 2 longitudal parameters ================
   static ROOT::Minuit2::Minuit2Minimizer minFitz(ROOT::Minuit2::kMigrad);

   Double_t R = TMath::Abs(trackPout[3]);
   for (int i = 0; i < RecoMiniMC_Nhits; i++) {
      Double_t X   = RecoMiniMC_XBHits[i];
      Double_t Y   = RecoMiniMC_YBHits[i];
      Double_t phi = TMath::ATan2(Y - yc, X - xc);
      if (RecoMiniMC_qHelix > 0)
         RecoMiniMC_tHelix[i] =
            (phi > phi0) ? twopi - RecoMiniMC_qHelix * (phi - phi0) : -RecoMiniMC_qHelix * (phi - phi0);
      else
         RecoMiniMC_tHelix[i] =
            (phi > phi0) ? -RecoMiniMC_qHelix * (phi - phi0) : twopi - RecoMiniMC_qHelix * (phi - phi0);
      /*printf("i=%2d R=(%6.1f %6.1f %6.1f) t=%6.1f phi0=%7.1f phi=%7.1f\n",i,
      RecoMiniMC_XBHits[i],RecoMiniMC_YBHits[i],RecoMiniMC_ZBHits[i],
      RecoMiniMC_tHelix[i],phi0*r2d,phi*r2d);*/
   }
   ROOT::Math::Functor fcn2(&RecoMiniMC_zHelix, 2);
   step[0] = 0.1;
   step[1] = 0.1;
   limU[0] = 100;
   limU[1] = 2000;
   limL[0] = -100;
   limL[1] = -2000;

   par[0] = 0; //
   par[1] = (RecoMiniMC_ZBHits[imax] - RecoMiniMC_ZBHits[imin]) /
            (RecoMiniMC_tHelix[imax] - RecoMiniMC_tHelix[imin]); // initial Z-velocity
   // printf("t1,t2,t(%6.1f %6.1f %6.1f) z2,z3(%6.1f %6.1f)\n",t1,t2,t,z2,z3);
   // printf("PinHelix: par(%8.1f %8.1f)\n",trackPin[2],trackPin[5]);
   // printf("fitHelix: par(%8.1f %8.1f)\n",par[0],par[1]);

   minFitz.SetFunction(fcn2);

   // Set the free variables to be minimized!
   minFitz.SetLimitedVariable(0, "x0", par[0], step[0], limL[0], limU[0]);
   minFitz.SetLimitedVariable(1, "y0", par[1], step[1], limL[1], limU[1]);

   res               = minFitz.Minimize();
   const double *xs2 = minFitz.X();
   if (minFitz.MinValue() / RecoMiniMC_Nhits > EvalpHit) return 2;
   // printf("outHelix: par(%8.1f %8.1f)\n",xs2[0],xs2[1]);
   // printf("chi2min2=%f\n",minFitz.MinValue()/RecoMiniMC_Nhits);
   //       minFitz.PrintResults();
   // cout<<"res="<<res<<endl;
   trackPout[2] = xs2[0];
   trackPout[5] = xs2[1];

   for (Int_t k = 0; k < 6; k++) trackDiff[k] = trackPout[k] - trackPin[k];

   //      pRecovery();
   return 0;
}

//______________________________________________________________________
void RecoMiniMC::CircleCenter(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t x3, Double_t y3,
                              Double_t &xc, Double_t &yc)
{
   Double_t a1, b1, c1, a2, b2, c2, d;

   a1 = x2 - x1;
   b1 = y2 - y1;
   a2 = x3 - x1;
   b2 = y3 - y1;
   c1 = 0.5 * (x2 * x2 - x1 * x1 + y2 * y2 - y1 * y1);
   c2 = 0.5 * (x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1);
   d  = a1 * b2 - a2 * b1;
   if (TMath::Abs(d) > 1e-12) { // crossing of 2 perpendiculars
      xc = (c1 * b2 - c2 * b1) / d;
      yc = (a1 * c2 - a2 * c1) / d;
   } else { // 3 points lie on a line
      TVector3 v  = TVector3(x3 - x1, y3 - y1, 0);
      TVector3 B  = TVector3(0, 0, 1);
      TVector3 rc = v.Cross(B);
      xc          = x1 + 100000 * rc.X() / rc.Mag();
      yc          = y1 + 100000 * rc.Y() / rc.Mag();
   }
}

//______________________________________________________________________
// put chi2 values
void RecoMiniMC::putChi2H()
{
   //    for(Int_t i=0;i<6;i++) trackPout[i]=trackPin[i];
   Double_t chisq = 0;
   Double_t delta;
   Double_t phi0 = trackPout[4];
   Double_t R    = TMath::Abs(trackPout[3]);
   Int_t    q    = trackPout[3] / R;
   for (int i = 0; i < RecoMiniMC_Nhits; i++) {
      Double_t t  = RecoMiniMC_tHelix[i];
      Double_t dx = RecoMiniMC_XBHits[i] - (trackPout[0] + R * TMath::Cos(phi0 - q * t));
      Double_t dy = RecoMiniMC_YBHits[i] - (trackPout[1] + R * TMath::Sin(phi0 - q * t));
      Double_t dz = RecoMiniMC_ZBHits[i] - (trackPout[2] + trackPout[5] * t);
      delta       = TMath::Sqrt(dx * dx + dy * dy + dz * dz) / RecoMiniMC_WHits[i];
      chisq += delta * delta;
      RecoMiniMC_chi2s[i] = delta * delta;
      fNhitSec[RecoMiniMC_iSec[i]]++;
      fSecChi2[RecoMiniMC_iSec[i]] += delta * delta;
   }
   Chi2Fit = chisq / RecoMiniMC_Nhits;
}

/*====================================================================
====================================================================
====================================================================
====================================================================


//______________________________________________________________________
  // Chi2 for the clusters fit by the helix
  Double_t RecoMiniMC_Chi2Helix(const Double_t *par){

    double chisq = 0;
    double t,delta;
    Int_t nhits=RecoMiniMC_Nhits;
const Double_t r2d=180/TMath::Pi();
//printf("nh=%d (%f,%f,%f,%f,%f,%f)\n",nhits,par[0],par[1],par[2],par[3],par[4],par[5]);
    double phi0=par[4];
    int q=RecoMiniMC_qHelix;
    for(int i=0;i<nhits;i++) {
      double X=RecoMiniMC_XBHits[i];
      double Y=RecoMiniMC_YBHits[i];
      double Z=RecoMiniMC_ZBHits[i];
      double W=RecoMiniMC_WHits[i];
      double dfi=phi0-TMath::ATan2(Y-par[1],X-par[0]);
      if(q>0) {
        if(dfi<0) t=TMath::TwoPi()+dfi; else t=dfi;
      }
      else {
        if(dfi>0) t=TMath::TwoPi()-dfi; else t=-dfi;
      }
      double dx=X-(par[0]+par[3]*TMath::Cos(par[4]-RecoMiniMC_qHelix*t));
      double dy=Y-(par[1]+par[3]*TMath::Sin(par[4]-RecoMiniMC_qHelix*t));
      double dz=Z-(par[2]+par[5]*t);
      delta=TMath::Sqrt(dx*dx+dy*dy+dz*dz)/W;
      chisq += delta*delta;
    }
//printf("   chi2=%f\n",chisq);
//exit(0);
    return chisq;
  }
//printf("nh=%d (%f,%f,%f,%f,%f,%f)\n",nhits,par[0],par[1],par[2],par[3],par[4]*r2,par[5]);
//printf("   chi2=%f\n",chisq);


//______________________________________________________________________
// Reconstruction of the helix parameters
  Int_t RecoMiniMC::fitHelix() {
  /*
    in the BCS the helix equation:
      x=p[0]+p[3]*cos(p[4]+q*t)
      y=p[1]+p[3]*sin(p[4]+q*t)
      z=p[2]+p[5]*t
  /
  // Convert global coordinates to MF coordinate system
  for(Int_t i=0;i<RecoMiniMC_Nhits;i++) {
    TVector3 rgClaster=TVector3(RecoMiniMC_XHits[i],RecoMiniMC_YHits[i],
                                RecoMiniMC_ZHits[i]);
    TVector3 rbClaster=G2B*rgClaster;
    RecoMiniMC_XBHits[i]=rbClaster.X();
    RecoMiniMC_YBHits[i]=rbClaster.Y();
    RecoMiniMC_ZBHits[i]=rbClaster.Z();
//printf("Ho(%f,%f,%f)\n", MpdTpcMiniMC_XHits[i],MpdTpcMiniMC_YHits[i],MpdTpcMiniMC_ZHits[i]);
//printf("H(%f,%f,%f)
HB(%f,%f,%f)\n",rgClaster.X(),rgClaster.Y(),rgClaster.Z(),rbClaster.X(),rbClaster.Y(),rbClaster.Z());
  }
    // initial approach
    Int_t imin,imax;
    Double_t rmin=200,rmax=0;
    for(Int_t i=0;i<MpdTpcMiniMC_Nhits;i++){
      Double_t rtpc=TMath::Sqrt(RecoMiniMC_XBHits[i]*RecoMiniMC_XBHits[i]+
                                RecoMiniMC_YBHits[i]*RecoMiniMC_YBHits[i]);
      if(rtpc<rmin) {imin=i;rmin=rtpc;}
      if(rtpc>rmax) {imax=i;rmax=rtpc;}
    }
    Double_t x1,y1,z1,a1,b1,c1,x2,y2,z2,a2,b2,c2,x3,y3,z3,
    xc,yc,xc1223,yc1223,xc1213,yc1213,xc1323,yc1323,d,t,t1,t2,t3,q;
    x2=RecoMiniMC_XBHits[imin];
    y2=RecoMiniMC_YBHits[imin];
    z2=RecoMiniMC_ZBHits[imin];
    x3=RecoMiniMC_XBHits[imax];
    y3=RecoMiniMC_YBHits[imax];
    z3=RecoMiniMC_ZBHits[imax];
    if(TMath::Abs(z3)<TMath::Abs(z2)) {
      x1=x2;y1=y2;z1=z2;
      x2=x3;y2=y3;z2=z3;
      x3=x2;y3=y2;z3=z2;
    }
    x1=0;
    y1=0;
    z1=0;

    a1=x2-x1;
    b1=y2-y1;
    a2=x3-x1;
    b2=y3-y1;
    c1=0.5*(x2*x2-x1*x1+y2*y2-y1*y1);
    c2=0.5*(x3*x3-x1*x1+y3*y3-y1*y1);
    d=a1*b2-a2*b1;
    if(TMath::Abs(d)>epsMin) {              // crossing of 2 perpendiculars
      xc=(c1*b2-c2*b1)/d;
      yc=(a1*c2-a2*c1)/d;
    }
    else {                                           // 3 points lie on a line
      TVector3 v=TVector3(x3-x1,y3-y1,0);
      TVector3 B=TVector3(0,0,1);
      TVector3 rc=v.Cross(B);
      xc=0.5*(x3+x1)+q*100000*rc.X()/rc.Mag();
      yc=0.5*(y3+y1)+q*100000*rc.Y()/rc.Mag();
    }

    t1=TMath::ATan2(y1-yc,x1-xc);
    t2=TMath::ATan2(y3-yc,x3-xc);
    if(t1>t2) {
      if(t1-t2>pi) {
        RecoMiniMC_qHelix=-1;
        t=twopi-(t1-t2);
      }
      else {
        RecoMiniMC_qHelix=1;
        t=t1-t2;
      }
    }
    else {
      if(t2-t1>pi) {
        RecoMiniMC_qHelix=1;
        t=twopi-(t2-t1);
      }
      else {
        RecoMiniMC_qHelix=-1;
        t=t2-t1;
      }
    }
    // Choose method upon creation between:
    // kMigrad, kSimplex, kCombined,
    // kScan, kFumili
    static ROOT::Minuit2::Minuit2Minimizer minFit(ROOT::Minuit2::kMigrad);

    ROOT::Math::Functor fcn(&RecoMiniMC_Chi2Helix,6);
    double step[6] = { 1,   1,  1,      1,         0.1,      1};
    double limU[6] = { 2000, 2000, 2000, 2000, TMath::Pi(), 2000};
    double limL[6] = {-2000,-2000,-2000,    0,-TMath::Pi(),-2000};
    double par[6];


    par[0]=xc;                  //
    par[1]=yc;                  //  - coordinates of the muon start point
    par[2]=0;                   //
    par[3]=TMath::Sqrt(xc*xc+yc*yc); // radius of the circle in XY plane of BCS
    par[4]=TMath::ATan2(-yc,-xc);       // initial muon angle (XY plane of BCS)
    par[5]=(z3-z1)/t;

/printf("fitHelix: par0(%8.1f %8.1f %4.1f %8.1f %6.1f %8.1f)  q=%d\n",
par[0],par[1],par[2],par[3],r2d*par[4],par[5],RecoMiniMC_qHelix);
printf("fitHelix: parIn(%8.1f %8.1f %4.1f %8.1f %6.1f %8.1f)\n",
trackPin[0],trackPin[1],trackPin[2],trackPin[3],r2d*trackPin[4],trackPin[5]);/

    minFit.SetFunction(fcn);

    // Set the free variables to be minimized!
    minFit.SetLimitedVariable(0,"x0",   par[0],step[0],limL[0],limU[0]);
    minFit.SetLimitedVariable(1,"y0",   par[1],step[1],limL[1],limU[1]);
    minFit.SetLimitedVariable(2,"z0",   par[2],step[2],limL[2],limU[2]);
    minFit.SetLimitedVariable(3,"R",    par[3],step[3],limL[3],limU[3]);
    minFit.SetLimitedVariable(4,"phi",  par[4],step[4],limL[4],limU[4]);
    minFit.SetLimitedVariable(5,"Pz",   par[5],step[5],limL[5],limU[5]);

    Bool_t res=minFit.Minimize();
//      minFit.PrintResults();
      const double *xs = minFit.X();
      if(minFit.MinValue()/RecoMiniMC_Nhits>EvalpHit) return 1;
      trackPout[0]=xs[0];
      trackPout[1]=xs[1];
      trackPout[2]=xs[2];
      trackPout[3]=RecoMiniMC_qHelix*xs[3];
      trackPout[4]=xs[4];
      trackPout[5]=xs[5];
//printf("fitHelix:parOut(%8.1f %8.1f %4.1f %8.1f %6.1f %8.1f)\n",
//trackPout[0],trackPout[1],trackPout[2],trackPout[3],r2d*trackPout[4],trackPout[5]);
      for(Int_t k=0;k<6;k++) trackDiff[k]=trackPout[k]-trackPin[k];
      return 0;

  }
//printf("fitHelix:  par(%8.1f %8.1f %4.1f %8.1f %6.1f %8.1f)\n",
//xs[0],xs[1],xs[2],xs[3],r2d*xs[4],xs[5]);
//printf("fitHelix:parIn(%8.1f %8.1f %4.1f %8.1f %6.1f %8.1f)\n",
//trackPin[0],trackPin[1],trackPin[2],trackPin[3],r2d*trackPin[4],trackPin[5]);

//______________________________________________________________________
  // put chi2 values
  void RecoMiniMC::putChi2H() {
    Double_t chisq = 0;
    Double_t delta;
    Int_t nhits=RecoMiniMC_Nhits;
    Double_t phi0=trackPout[4];
    Double_t R=TMath::Abs(trackPout[3]);
    Int_t q=trackPout[3]/R;
//printf("putChi2: trackPin=(%f,%f,%f,%f,%f,%f,)\n",
//trackPin[0],trackPin[1],trackPin[2],trackPin[3],trackPin[4],trackPin[5]);
//printf("putChi2: trackPout=(%f,%f,%f,%f,%f,%f,)\n",
//trackPout[0],trackPout[1],trackPout[2],trackPout[3],trackPout[4],trackPout[5]);
    for(int i=0;i<nhits;i++) {
      Double_t X=RecoMiniMC_XBHits[i];
      Double_t Y=RecoMiniMC_YBHits[i];
      Double_t Z=RecoMiniMC_ZBHits[i];
      Double_t W=RecoMiniMC_WHits[i];
      Double_t phi=TMath::ATan2(Y-trackPout[1],X-trackPout[0]);
      Double_t t;
      if(q>0) t=(phi>phi0)?twopi-q*(phi-phi0):-q*(phi-phi0);
      else    t=(phi>phi0)?-q*(phi-phi0):twopi-q*(phi-phi0);
      Double_t dx=X-(trackPout[0]+R*TMath::Cos(phi0-q*t));
      Double_t dy=Y-(trackPout[1]+R*TMath::Sin(phi0-q*t));
      Double_t dz=Z-(trackPout[2]+trackPout[5]*t);
      delta=TMath::Sqrt(dx*dx+dy*dy+dz*dz)/W;
      chisq += delta*delta;
      RecoMiniMC_chi2s[i]= delta*delta;
      chisq += delta*delta;
      fNhitSec[RecoMiniMC_iSec[i]]++;
      fSecChi2[RecoMiniMC_iSec[i]]+=delta*delta;
    }
    Chi2Fit=chisq/RecoMiniMC_Nhits;
  }
//      if(delta*delta>1000)
//printf("putChi2: RecoMiniMC_iSec[%d]=%d delta=%f  chi2=%f\n",i,RecoMiniMC_iSec[i],delta*delta,chisq);
//printf("putChi2: RecoMiniMC_iSec[%d] delta^2=%f  chi2=%f\n",i,delta*delta,chisq);
//      if(Chi2Fit>100)
//printf("putChi2: hits=%d Chi2Fit per hit = %f\n",RecoMiniMC_Nhits,Chi2Fit);
====================================================================
====================================================================
====================================================================
====================================================================*/

//______________________________________________________________________

// create new cluster
void RecoMiniMC::convert()
{

   /*
      <Xgo,Ygo,Zgo> - global coordinates of the hit(old alignment)
      <xlo,ylo,zlo> - local coordinates of the hit(old alignment)
      <Xl,Yl> - local coordinates of the cluster
      <xln,Yln,zln> - local coordinates of the hit(new alignment)
      <Xgn,Ygn,Zgn> - global coordinates of the hit(new alignment)
      <xlo,ylo,zlo>=<G0_gsc>+{glo2loc}*<Xgo,Ygo,Zgo>
      <e_3o> - unit vector in GSC along the Z-axis of the old LSC
      <e_3n> - unit vector in GSC along the Z-axis of the new LSC
      <eef> - unit vector in GSC of the electric field
      {glo2loc} - transformation matrix GSC->LSC
      {loc2glo} - transformation matrix LSC->GSC
      cosAB=(e_3o*eef)  Zalong=-zlo/cosAB <Bo>=-zlo*<e_3o> <A>=Zalong*eef
      <AmB>=<A>-<BO>  <amb>={glo2loc}_old*<AmB>
         Xl=xlo+<amb>_x  Yl=xlo+<amb>_y
       - - - - - - new alignment - - - -
      cosAB=(e_3n*eef) <Bn>=(-Zalong*cosAB)*<e_3n>
      <AmB>=<A>-<Bn>  <amb>={glo2loc}_new*<AmB>
         xln=Xl-<amb>_x  xln=Xl-<amb>_y  zln=<Bn>_x
         <Xgn,Ygn,Zgn>=R0shift+{loc2glo}_new*<xln,Yln,zln>
   */
   Double_t cosAB, ZalongE, XlClaster, YlClaster;
   TVector3 eef, e_3n, e_3o, A, B, AmB, amb, rghit, rlhit, rbhit;
   Int_t    is, iso = 100;
   RecoMiniMC_Nhits = 0;
   for (Int_t i = 0; i < MpdTpcMiniMC_Nhits; i++) {
      is = MpdTpcMiniMC_iSec[i];
      if (sectors[is]) {
         if (is != iso) {
            if (is < 12) {
               if (eEf.Z() > 0)
                  eef = eEf.Unit();
               else
                  eef = -eEf.Unit();
            } else {
               if (eEf.Z() > 0)
                  eef = -eEf.Unit();
               else
                  eef = eEf.Unit();
            }
            R0shift  = newR0_shift[is];
            loc2glo  = new_loc2glo[is];
            G0_gsc   = newG0_gsc[is];
            glo2loc  = new_glo2loc[is];
            oG0_gsc  = oldG0_gsc[is];
            oglo2loc = old_glo2loc[is];
            e_3o     = TVector3(oglo2loc.ZX(), oglo2loc.ZY(), oglo2loc.ZZ());
            e_3n     = TVector3(loc2glo.XZ(), loc2glo.YZ(), loc2glo.ZZ());
            // printf("e_3o(%9.4f,%9.4f,%9.4f) e_3n(%9.4f,%9.4f,%9.4f) sec=%d\n",
            // e_3o.X(),e_3o.Y(),e_3o.Z(),e_3n.X(),e_3n.Y(),e_3n.Z(),is);
            iso = is;
         }
         rghit       = TVector3(MpdTpcMiniMC_XHits[i], MpdTpcMiniMC_YHits[i],
                                MpdTpcMiniMC_ZHits[i]); // hit in GCS
         TVector3 ho = rghit;
         rlhit       = oG0_gsc + oglo2loc * rghit; // hit in LCS
         cosAB       = TMath::Abs(e_3o * eef);
         ZalongE     = TMath::Abs(rlhit.Z()) / cosAB; // distance: track->cluster
         A           = ZalongE * eef;                 // vector: track->cluster
         B           = (ZalongE * cosAB) * e_3o;      // perpendicular from trak point to the cluster
         AmB         = A - B;                         // vector in the sector plane: perpendicular base->cluster(GCS)
         amb         = oglo2loc * AmB;                // vector AmB in LCS
         XlClaster   = rlhit.X() + amb.X();
         YlClaster   = rlhit.Y() + amb.Y();
         //  - - - - - - new alignment - - - -
         cosAB                              = e_3n * eef;
         B                                  = (ZalongE * cosAB) * e_3n;
         AmB                                = A - B;
         amb                                = glo2loc * AmB;
         rlhit                              = TVector3(XlClaster - amb.X(), YlClaster - amb.Y(), -ZalongE * cosAB);
         rghit                              = R0shift + loc2glo * rlhit;
         RecoMiniMC_XHits[RecoMiniMC_Nhits] = rghit.X();
         RecoMiniMC_YHits[RecoMiniMC_Nhits] = rghit.Y();
         RecoMiniMC_ZHits[RecoMiniMC_Nhits] = rghit.Z();
         if (MFieldMode) {
            rbhit                               = G2B * rghit; // GCS->BCS
            RecoMiniMC_XBHits[RecoMiniMC_Nhits] = rbhit.X();
            RecoMiniMC_YBHits[RecoMiniMC_Nhits] = rbhit.Y();
            RecoMiniMC_ZBHits[RecoMiniMC_Nhits] = rbhit.Z();
         }
         RecoMiniMC_WHits[RecoMiniMC_Nhits] = MpdTpcMiniMC_WHits[i];
         RecoMiniMC_iSec[RecoMiniMC_Nhits]  = is;
         RecoMiniMC_Nhits++;
         // printf("h=%d rIn(%9.4f,%9.4f,%9.4f) rOut(%9.4f,%9.4f,%9.4f) sec=%d\n",
         // i,rlhit.X(),rlhit.Y(),rlhit.Z(),rghit.X(),rghit.Y(),rghit.Z(),is);
         // printf("     oh(%9.4f,%9.4f,%9.4f) nh(%9.4f,%9.4f,%9.4f) sec=%d\n",
         // ho.X(),ho.Y(),ho.Z(),rghit.X(),rghit.Y(),rghit.Z(),is);

         // printf("i=%2d old(%9.4f,%9.4f,%9.4f) new(%9.4f,%9.4f,%9.4f)\n",i,
         // MpdTpcMiniMC_XHits[i],MpdTpcMiniMC_YHits[i],MpdTpcMiniMC_ZHits[i],
         // RecoMiniMC_XHits[i],RecoMiniMC_YHits[i],RecoMiniMC_ZHits[i]);
      }
   }
}

//______________________________________________________________________
