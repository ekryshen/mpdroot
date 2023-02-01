///
/// \brief MiniMC for TPC alignment tests
///
/// \author Valentin Kuzmin, SINP of Moscow State University

#include "TpcMiniMC.h"
#include <iostream>
#include <cstdlib>
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
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

// Collaborating Class Headers --------

//_____________________________________________________________________

// Constructor 2 difines all parameters
TpcMiniMC::TpcMiniMC(BaseTpcSectorGeo &tpcGeo, Double_t *B_mf, Bool_t Beam_point, Double_t *E_ef, Double_t *P_las,
                     Double_t D_track, Double_t *Min_padS, Double_t *Max_padS, Double_t *Sig_padS,
                     Double_t Sig_timeBins, Double_t Sig_Z_pos, Double_t *R0_beam, Double_t K_par, Double_t P_min,
                     Double_t P_max, Int_t Chambers, Int_t MaxPcluster, Int_t MinHits, UInt_t Seed, TString *aFileIn,
                     TString *OutDir, TString *OutFile = new TString(" "))
{
   // true if the Magnetic field is on
   MFieldMode = TMath::Abs(B_mf[0]) > epsMin || TMath::Abs(B_mf[1]) > epsMin || TMath::Abs(B_mf[2]) > epsMin;
   // store parameters
   chambers = Chambers;
   if (Chambers == 0) {
      sector0 = 0;
      sector1 = 24;
   } else if (Chambers == 1) {
      sector0 = 0;
      sector1 = 12;
   } else if (Chambers == 2) {
      sector0 = 0;
      sector1 = 24;
   } else
      cout << " MiniMC: no chambers mode " << Chambers << endl;
   fminHits   = MinHits;
   beam_point = Beam_point;
   for (Int_t i = 0; i < 3; i++) {
      e_ef[i]    = E_ef[i];
      b_mf[i]    = B_mf[i];
      r0_beam[i] = R0_beam[i];
      p_las[i]   = P_las[i];
   }
   for (Int_t i = 0; i < 2; i++) {
      min_padS[i] = Min_padS[i];
      max_padS[i] = Max_padS[i];
      sig_padS[i] = Sig_padS[i];
   }
   sig_timeBins = Sig_timeBins;
   sig_Z_pos    = Sig_Z_pos;
   k_par        = K_par;
   // printf("TpcMiniMC: k_par=%d\n",k_par);
   p_min       = P_min;
   p_max       = P_max;
   maxPcluster = MaxPcluster;
   d_track     = D_track;
   fSeed       = Seed;
   aFile       = aFileIn;
   outDir      = OutDir;
   outFile     = OutFile;

   Init(tpcGeo);
}

// Constructor 3 difines reads parameters from the file
TpcMiniMC::TpcMiniMC(BaseTpcSectorGeo &tpcGeo, UInt_t Seed, TString *ParFile, TString *aFileIn, TString *OutDir,
                     TString *OutFile)
{
   aFile   = aFileIn;
   outDir  = OutDir;
   outFile = OutFile;
   parFile = ParFile;

   const char *pFile   = parFile->Data();
   TFile      *parF    = new TFile(pFile, "READ");
   TTree      *t_parIn = (TTree *)parF->Get("parameters");

   t_parIn->SetBranchAddress("b_mf", &b_mf);
   t_parIn->SetBranchAddress("e_ef", &e_ef);
   t_parIn->SetBranchAddress("r0_beam", &r0_beam);
   t_parIn->SetBranchAddress("p_las", &p_las);
   t_parIn->SetBranchAddress("d_track", &d_track);
   t_parIn->SetBranchAddress("p_min", &p_min);
   t_parIn->SetBranchAddress("p_max", &p_max);
   t_parIn->SetBranchAddress("min_padS", &min_padS);
   t_parIn->SetBranchAddress("max_padS", &max_padS);
   t_parIn->SetBranchAddress("sig_padS", &sig_padS);
   t_parIn->SetBranchAddress("sig_timeBins", &sig_timeBins);
   t_parIn->SetBranchAddress("sig_Z_pos", &sig_Z_pos);
   t_parIn->SetBranchAddress("sector0", &sector0);
   t_parIn->SetBranchAddress("sector1", &sector1);
   t_parIn->SetBranchAddress("k_par", &k_par);
   t_parIn->SetBranchAddress("maxPcluster", &maxPcluster);
   t_parIn->SetBranchAddress("fSeed", &fSeed);
   t_parIn->SetBranchAddress("beam_point", &beam_point);
   t_parIn->SetBranchAddress("MFieldMode", &MFieldMode);
   t_parIn->GetEvent(0);
   parF->Close();
   fSeed = Seed;
   if (sector0 == 0) {
      if (sector1 == 12)
         chambers = 1;
      else
         chambers = 0;
   } else
      chambers = 2;

   Init(tpcGeo);
}

//_____________________________________________________________________
void TpcMiniMC::Init(BaseTpcSectorGeo &tpcGeo)
{
   //
   if (chambers == 0) {
      sector0 = 0;
      sector1 = 24;
   } else {
      if (chambers == 1) {
         sector0 = 0;
         sector1 = 12;
         j1_ls   = 0;
      } else {
         if (chambers == 2) {
            sector0 = 12;
            sector1 = 24;
            j2_ls   = 0;
         } else {
            cout << " MiniMC: no chambers mode " << chambers << endl;
            return;
         }
      }
   }
   Double_t mod_p = sqrt(p_las[0] * p_las[0] + p_las[1] * p_las[1] + p_las[2] * p_las[2]);
   for (Int_t i = 0; i < 3; i++) p_las[i] = p_las[i] / mod_p;
   // - - - - - - -

   fTpcSecGeo = dynamic_cast<TpcSectorGeoAlignmentVAK *>(&tpcGeo);
   if (!fTpcSecGeo) Fatal("MpdTpcMiniMC::Init", " !!! Wrong geometry type !!! ");
   // - - - - - - -
   if (!aFile->EqualTo(" ")) { // read a special alignment
      const char *aInFile = aFile->Data();
      alignmentF          = new TFile(aInFile, "READ");
      t_alignment         = (TTree *)alignmentF->Get("alignment");
      t_alignment->SetBranchAddress("R0_A", &R0_A);
      t_alignment->SetBranchAddress("alpha_A", &alpha_A);
      t_alignment->SetBranchAddress("beta_A", &beta_A);
      t_alignment->SetBranchAddress("gamma_A", &gamma_A);
      t_alignment->GetEvent(0);
      alignmentF->Close();
      fTpcSecGeo->TpcSectorGeoA(R0_A, alpha_A, beta_A, gamma_A);
   }
   for (Int_t i = 0; i < 24; i++) {
      TVector3 r0l;
      fTpcSecGeo->GetAlignment(i, r0l, alpha_A[i], beta_A[i], gamma_A[i]);
      R0_A[0][i] = r0l.X();
      R0_A[1][i] = r0l.Y();
      R0_A[2][i] = r0l.Z();
   }
   t_alignment = (TTree *)alignmentF->Get("alignment");
   fL_svol     = fTpcSecGeo->zPads().X();
   fZmin       = fTpcSecGeo->GetZmin();
   fZmax       = fTpcSecGeo->GetZmax();

   sZ_min = fZmax; // using in Helix z calculation

   bMf     = TVector3(b_mf[0], b_mf[1], b_mf[2]);
   bMF_mod = bMf.Mag();
   eEf     = TVector3(e_ef[0], e_ef[1], e_ef[2]);
   r0Beam  = TVector3(r0_beam[0], r0_beam[1], r0_beam[2]);
   pLas    = TVector3(p_las[0], p_las[1], p_las[2]);

   // set magnetic field parameters
   eBz              = bMf.Unit();
   Double_t eGyeBz  = eGy.Dot(eBz);
   Double_t weGyeBz = 1 / TMath::Sqrt(1 - eGyeBz);
   eBy              = weGyeBz * (eGy + eGyeBz * eBz);
   eBx              = eBy.Cross(eBz);
   Double_t phi_B;
   if (bMf.X() * bMf.X() + bMf.Y() * bMf.Y() == 0)
      phi_B = 0;
   else
      phi_B = TMath::ATan2(bMf.Y(), bMf.X());
   B2G.SetToIdentity();
   B2G.RotateZ(phi_B);
   Double_t theta_B = TMath::ACos(bMf.Z() / bMF_mod);
   B2G.RotateY(theta_B);
   G2B         = B2G.Inverse();
   Int_t nRows = fTpcSecGeo->NofRows();
   // fill min&max pads coordinates
   for (Int_t irow = 0; irow < nRows; irow++) {
      Int_t nPads = fTpcSecGeo->NofPadsInRow(irow);
      for (Int_t ipad = 0; ipad < nPads; ipad++) {
         fTpcSecGeo->MinMaxPadPosition(irow, ipad, xminPad[ipad][irow], yminPad[ipad][irow], xmaxPad[ipad][irow],
                                       ymaxPad[ipad][irow]);
      }
   }

   // set names of header&data files
   stringstream ss;
   //  time_t now = time(NULL);
   //  ss<<"miniMC"<<now<<".root";
   if (!outFile->EqualTo(" ")) { // read a special alignment
      ss << outDir->Data() << "/" << outFile->Data() << ".root";
   } else {
      Ssiz_t  ld     = aFile->Last('/');
      Ssiz_t  ln     = aFile->Last('.');
      TString fName0 = aFile->Data();
      TString fName  = fName0(ld + 1, ln - ld - 1);
      ss << outDir->Data() << "/miniMC_" << fName.Data() << ".root";
   }
   cout << "miniMC output file name: " << ss.str() << endl;
   TString *oFile = new TString(ss.str());
   outfile        = oFile->Data();

   outmMCfile = new TFile(outfile, "RECREATE");
   fRandom->SetSeed(fSeed);

   t_parameters = new TTree("parameters", "miniMC_parameters");

   t_parameters->Branch("b_mf", &b_mf, "b_mf[3]/D");
   t_parameters->Branch("e_ef", &e_ef, "e_ef[3]/D");
   t_parameters->Branch("r0_beam", &r0_beam, "r0_beam[3]/D");
   t_parameters->Branch("p_las", &p_las, "p_las[3]/D");
   t_parameters->Branch("d_track", &d_track, "d_track/D");
   t_parameters->Branch("p_min", &p_min, "p_min/D");
   t_parameters->Branch("p_max", &p_max, "p_max/D");
   t_parameters->Branch("min_padS", &min_padS, "min_padS[2]/D");
   t_parameters->Branch("max_padS", &max_padS, "max_padS[2]/D");
   t_parameters->Branch("sig_padS", &sig_padS, "sig_padS[2]/D");
   t_parameters->Branch("sig_timeBins", &sig_timeBins, "sig_timeBins/D");
   t_parameters->Branch("sig_Z_pos", &sig_Z_pos, "sig_Z_pos/D");
   t_parameters->Branch("sector0", &sector0, "sector0/I");
   t_parameters->Branch("sector1", &sector1, "sector1/I");
   t_parameters->Branch("k_par", &k_par, "k_par/I");
   t_parameters->Branch("maxPcluster", &maxPcluster, "maxPcluster/I");
   t_parameters->Branch("fSeed", &fSeed, "fSeed/i");
   t_parameters->Branch("beam_point", &beam_point, "beam_point/O");
   t_parameters->Branch("MFieldMode", &MFieldMode, "MFieldMode/O");

   t_alignment = new TTree("alignment", "alignment_parameters");
   t_alignment->Branch("R0_A", &R0_A, "R0_A[3][24]/D");
   t_alignment->Branch("alpha_A", &alpha_A, "alpha_A[24]/D");
   t_alignment->Branch("beta_A", &beta_A, "beta_A[24]/D");
   t_alignment->Branch("gamma_A", &gamma_A, "gamma_A[24]/D");
   t_alignment->Fill();
   t_alignment->Write();

   // Get global sctor parameters
   fWpad        = fTpcSecGeo->PadWidth(0);
   fHpad[0]     = fTpcSecGeo->PadHeight(0);
   fHpad[1]     = fTpcSecGeo->PadHeight(27);
   fNrows[0]    = fTpcSecGeo->NofRowsReg(0);
   fNrows[1]    = fTpcSecGeo->NofRowsReg(1);
   fNrows[2]    = fTpcSecGeo->NofRows();
   dh_track     = d_track / 2;
   minMuonTheta = TMath::ATan(45. / 170.);

   t_data = new TTree("data", "miniMC_data");
   t_data->Branch("trackPin", &trackPin, "trackPin[6]/D");
   t_data->Branch("trackPout", &trackPout, "trackPout[6]/D");
   t_data->Branch("Chi2Fit", &Chi2Fit, "Chi2Fit/D");
   t_data->Branch("Par0R", &trackDiff[0], "trackDiff/D");
   t_data->Branch("Par1R", &trackDiff[1], "trackDiff/D");
   t_data->Branch("Par2R", &trackDiff[2], "trackDiff/D");
   t_data->Branch("Par3R", &trackDiff[3], "trackDiff/D");
   t_data->Branch("Par4R", &trackDiff[4], "trackDiff/D");
   t_data->Branch("Par5R", &trackDiff[5], "trackDiff/D");
   t_data->Branch("MpdTpcMiniMC_Nhits", &MpdTpcMiniMC_Nhits, "MpdTpcMiniMC_Nhits/I");
   t_data->Branch("MpdTpcMiniMC_iSec", &MpdTpcMiniMC_iSec, "MpdTpcMiniMC_iSec[MpdTpcMiniMC_Nhits]/I");
   t_data->Branch("MpdTpcMiniMC_XHits", &MpdTpcMiniMC_XHits, "MpdTpcMiniMC_XHits[MpdTpcMiniMC_Nhits]/D");
   t_data->Branch("MpdTpcMiniMC_YHits", &MpdTpcMiniMC_YHits, "MpdTpcMiniMC_YHits[MpdTpcMiniMC_Nhits]/D");
   t_data->Branch("MpdTpcMiniMC_ZHits", &MpdTpcMiniMC_ZHits, "MpdTpcMiniMC_ZHits[MpdTpcMiniMC_Nhits]/D");
   t_data->Branch("MpdTpcMiniMC_WHits", &MpdTpcMiniMC_WHits, "MpdTpcMiniMC_WHits[MpdTpcMiniMC_Nhits]/D");
   t_data->Branch("MpdTpcMiniMC_chi2s", &MpdTpcMiniMC_chi2s, "MpdTpcMiniMC_chi2s[MpdTpcMiniMC_Nhits]/D");

   // fullChi2 TRee
   tout_fullChi2 = new TTree("Chi2", "fullChi2");

   tout_fullChi2->Branch("fNhits", &fNhits, "fNhits/L");
   tout_fullChi2->Branch("fSumChi2", &fSumChi2, "fSumChi2/D");
   tout_fullChi2->Branch("secNhits", &fNhitSec, "fNhitSec[24]/L");
   tout_fullChi2->Branch("secSumChi2", &fSecChi2, "fSecChi2[24]/D");

   ls1   = true;
   ls2   = false;
   i_ls1 = i1_ls;
   j_ls1 = j1_ls, a_ls1 = a11_ls;

   // cout<<"MpdTpcminiMC: Initialization  finished"<<endl;
}

//______________________________________________________________________
// simulate fnEvents and write data to OutPut file
void TpcMiniMC::simulate(Int_t nEvents)
{
   if (printL > 0) {
      if (!aFile->EqualTo(" ")) // read a special alignment
         cout << "MpdTpcMiniMC::simulate: Alignment data will be taken "
              << "from the file: \n"
              << aFile->Data() << endl;
      else
         cout << "MpdTpcMiniMC::simulate: Zero alignment data will be used " << endl;
      cout << "MpdTpcMiniMC::simulate: " << nEvents << " are to simulate\n";
   }
   fNevents       = nEvents;
   Int_t iEvt     = 0;
   Nreconstructed = 0;
   for (Int_t is = 0; is < 24; is++) {
      fNhitSec[is] = 0;
      fSecChi2[is] = 0;
   }
   while (Nreconstructed < fNevents) {
      // printf("simulate:  event=%d =======================\n",iEvt);
      iEvt++;
      Int_t attempts = miniEvent();
      if (attempts > 0 && MpdTpcMiniMC_Nhits >= fminHits) {
         if (printL > 1)
            if (iEvt % printL == 0)
               cout << "MpdTpcMiniMC: event=" << iEvt << " simulated on the " << attempts << " (th) attempt" << endl;
         if (!MFieldMode) {
            if (fitLine() != 0) continue;
            putChi2L();
         } else {
            // printf("   simulate: call fitHelix(%d) =============\n",iEvt);
            if (fitHelix() != 0) continue;
            putChi2H();
         }
         Nreconstructed++;
         // printf("   simulate: call putChi2(%d) =============\n",MpdTpcMiniMC_Nhits);
         t_data->Fill();
      } else {
         if (attempts < 0) {
            cout << "MpdTpcMiniMC: event=" << iEvt + 1 << " was not simulated on the " << attempts << " (th) attempt"
                 << endl;
            return;
         }
      }
      if (iEvt > fNevents && iEvt > 3 * Nreconstructed) {
         printf("STOP: very low efficienty of %6.2f at %d simulated events\n",
                float_t(Nreconstructed) / float_t(iEvt) * 100, iEvt);
         return;
      }
   }
   t_data->Write();
   printf("******  %d events reconstucted of %d simulated (%%%4.1f) ******\n", Nreconstructed, iEvt,
          float_t(Nreconstructed) / float_t(iEvt) * 100);
   fSeed = fRandom->GetSeed();
   t_parameters->Fill();
   t_parameters->Write();

   for (Int_t is = 0; is < 24; is++) {
      fNhits += fNhitSec[is];
      fSumChi2 += fSecChi2[is];
   }
   if (fNhits > 0) fSumChi2 /= fNhits;
   tout_fullChi2->Fill();
   tout_fullChi2->Write();

   outmMCfile->Close();
   if (printL > 0) cout << " =============  FinalSeed = " << fSeed << endl;
}

//______________________________________________________________________
// an event generation
Int_t TpcMiniMC::miniEvent()
{
   // printf("===== event=%d ================================\n",Nreconstructed);
   for (Int_t n_sim = 1; n_sim <= maxSim; n_sim++) {
      if (!MFieldMode) {
         if (beam_point) {
            if (k_par == 3) { // TPC laser system
               GetLaserRay();
            } else {
               if (k_par == 4) { // random muon
                  GetBeamMuon();
               } else {
                  if (k_par != 0) { // one direction partile from a fixed point
                     cout << "MpdTpcMiniMC::miniEvent: No simulation mode with" << endl;
                     cout << "  MFieldMode=" << MFieldMode << "  beam_point=" << beam_point << "  k_par=" << k_par
                          << endl;
                     return -2;
                  }
               }
            }
         } else {
            GetCosmicMuon();
         }
         crossLinePoints();
         // printf("miniEvent: global track (%f,%f,%f) -> (%f,%f,%f)\n",
         // fP11.X(),fP11.Y(),fP11.Z(),fP12.X(),fP12.Y(),fP12.Z());
         if (laser_ray() > 0) {
            return n_sim;
         }
      } else {
         if (beam_point) { // particles produce in one point
            GetPointMuon();
            if (passMuoninMF() > 0) {
               return n_sim;
            }
         } else {
            cout << "MpdTpcMiniMC::miniEvent: No simulation mode with" << endl;
            cout << "  MFieldMode=" << MFieldMode << "  beam_point=" << beam_point << "  k_par=" << k_par << endl;
            return -2;
         }
      }
   }
   return -1; // number simulation attemps was exeeded
}

//______________________________________________________________________

// laser ray event simulation
Int_t TpcMiniMC::laser_ray()
{
   MpdTpcMiniMC_Nhits = 0;
   for (Int_t isec = sector0; isec < sector1; isec++) {
      if ((isec < 12 && fP11.Mag2() == 0) || (isec > 11 && fP21.Mag2() == 0))
         continue;
      else { // There is an track in the sector sensitive region
         iSector = isec;
         if (isec < 12) {
            if (!crossCheck(isec, fP11, fP12)) continue;
         } else {
            if (!crossCheck(isec, fP21, fP22)) continue;
         }
      }
      // printf("laser_ray: sector=%d local track (%f,%f) -> (%f,%f)\n",
      // isec,fR0sec.X(),fR0sec.Y(),fR0sec1.X(),fR0sec1.Y());
      lineClusters(isec);
   }
   // printf("laser_ray: MpdTpcMiniMC_Nhits=%d\n",MpdTpcMiniMC_Nhits);
   return MpdTpcMiniMC_Nhits;
}

//______________________________________________________________________
// find clusters of direct line track
Int_t TpcMiniMC::lineClusters(Int_t isec)
{

   Int_t    iRow, iPad, nPads, iPad1, iPad2, nclust;
   Int_t    iRow1, iRow0, iRow2;
   Double_t xloc1, yloc1, xloc2, yloc2, yloc0, ylmin, ylmax, XPads, x1pad, x2pad;
   Double_t UpRowPart, LowRowPart, FirstPadPart, LastPadPart, RowPart, PadPart, spad;
   Double_t hPad, wPad;
   wPad   = fWpad;
   nclust = 0; // number of found clusters
   xloc1  = fR0sec.X();
   if (xloc1 > fR0sec1.X()) {
      xloc2 = xloc1;
      xloc1 = fR0sec1.X();
   } else {
      xloc2 = fR0sec1.X();
   }
   // if(Nreconstructed==11&&isec==2)
   // printf("lineClusters: fR1sec.y(%f->%f) fR2sec.y(%f->%f) tng=%f\n",
   // fR1sec.Y(),fR1sec1.Y(),fR2sec.Y(),fR2sec1.Y(),fEsec.X()/fEsec.Y());
   yloc1 = fR1sec.Y() < fR2sec.Y() ? fR1sec.Y() : fR2sec.Y();
   yloc2 = fR1sec1.Y() > fR2sec1.Y() ? fR1sec1.Y() : fR2sec1.Y();
   if (yloc1 > yloc2) {
      Double_t w = yloc1;
      yloc1      = yloc2;
      yloc2      = w;
   }
   if (yloc1 > fTpcSecGeo->GetMaxY() || yloc2 < 0) {
      return 0; // the track doesn't cross the sector
   }
   XPads = -fTpcSecGeo->GetMaxX();
   if (xloc1 > -XPads || xloc2 < XPads) {
      return 0; // the track doesn't cross the sector
   }
   // find the parts of bottom and top rows not covered by the track band
   iRow1 = fTpcSecGeo->NofRow(yloc1);
   iRow2 = fTpcSecGeo->NofRow(yloc2);
   if (iRow1 < 0) iRow1 = 0;
   if (iRow2 < 0) iRow2 = fNrows[2] - 1;
   // printf("lineClusters: s=%d iRow1,2(%d %d) xloc1,2(%f,%f)  yloc1,2(%f,%f)\n",
   // isec,iRow1,iRow2,xloc1,xloc2,yloc1,yloc2);
   ylmax = fTpcSecGeo->UpRowEdge(iRow2);
   hPad  = fTpcSecGeo->PadHeight(iRow2);
   if (yloc2 >= ylmax) {
      UpRowPart = 1;
   } else {
      UpRowPart = (yloc2 - ylmax + hPad) / hPad;
   }
   ylmin = fTpcSecGeo->LowRowEdge(iRow1);
   hPad  = fTpcSecGeo->PadHeight(iRow1);
   if (yloc1 <= ylmin) {
      LowRowPart = 1;
   } else {
      LowRowPart = (ylmin + hPad - yloc1) / hPad;
   }
   // printf("lineClusters: 1)hPad=%f ylmin=%f ylmax=%f LowPart=%f UpPart=%f \n",
   // hPad,ylmin,ylmax,LowRowPart,UpRowPart);
   //  - - - - - - - - - - - - - -
   //  - - - - - - - - - - - - - -
   if (TMath::Abs(fEsec.Y()) < epsMin) { //************ the track is || to X-axes
      // cicle over all rows crossed by the track
      for (Int_t irow = iRow1; irow <= iRow2; irow++) {
         nPads = fTpcSecGeo->NofPadsInRow(irow);
         XPads = -nPads * fWpad / 2;
         // printf("lineClusters: XPads=%f\n",XPads);
         //  find the parts of 1st and tlast pad not covered by the track band
         if (xloc1 > -XPads || xloc2 < XPads) {
            return 0; // the track doesn't cross the row
         }
         if (xloc1 <= XPads) { // the track starts befor the 1st pad
            FirstPadPart = 1.;
            iPad1        = 0;
         } else { // the track starts inside the row
            iPad1        = Int_t((xloc1 - XPads) / fWpad);
            FirstPadPart = xloc1 - (XPads + fWpad * iPad1);
         }
         if (xloc2 >= -XPads) { // the track stops behind the last pad
            LastPadPart = 1.;
            iPad2       = nPads - 1;
         } else { // the track stops inside the row
            iPad2       = Int_t((xloc2 - XPads) / fWpad);
            LastPadPart = XPads + fWpad * (iPad2 + 1) - xloc2;
         }
         //  find the row height covered by the track
         hPad = fTpcSecGeo->PadHeight(irow);
         if (irow == iRow1 && irow == iRow2)
            RowPart = LowRowPart + UpRowPart - 1;
         else if (irow == iRow1)
            RowPart = LowRowPart;
         else if (irow == iRow2)
            RowPart = UpRowPart;
         else
            RowPart = 1;
         // printf("lineClusters: 2)hPad=%f RowPart=%f\n",hPad,RowPart);
         for (Int_t ipad = iPad1; ipad <= iPad2; ipad++) {
            //  check cover of 1st and last pads
            if (ipad == iPad1 && ipad == iPad2)
               PadPart = FirstPadPart + LastPadPart - 1;
            else if (ipad == iPad1)
               PadPart = FirstPadPart;
            else if (ipad == iPad2)
               PadPart = LastPadPart;
            else
               RowPart = 1;
            SPads[ipad] = wPad * d_track * RowPart * PadPart;
            // printf("lineClusters: %d-%d-%d   S=%f (%f,%f)\n",isec,irow,ipad,SPads[ipad],wPad,hPad);
         }
         Int_t iRoc = fTpcSecGeo->NofRoc(irow);
         nclust += ClustersInRow(irow, iPad1, iPad2, iRoc, SPads);
      }
      return nclust;
   } // end  of the case: the track is || to the sector base

   // - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - -

   if (TMath::Abs(fEsec.X()) < epsMin) { //************ track is || to the y-axes
      xloc1 = fR1sec.X();
      if (xloc1 < fR2sec.X()) { // find min&max X of the band
         xloc2 = fR2sec.X();
      } else {
         xloc2 = xloc1;
         xloc1 = fR2sec.X();
      }
      for (Int_t irow = iRow1; irow <= iRow2; irow++) {
         nPads      = fTpcSecGeo->NofPadsInRow(irow);
         XPads      = -nPads * fWpad / 2;
         Int_t iRoc = fTpcSecGeo->NofRoc(irow);
         if (xloc2 < XPads || xloc1 > -XPads) continue;
         hPad = fTpcSecGeo->PadHeight(irow);
         if (irow == iRow1 && irow == iRow2)
            RowPart = LowRowPart + UpRowPart - 1;
         else if (irow == iRow1)
            RowPart = LowRowPart;
         else if (irow == iRow2)
            RowPart = UpRowPart;
         else
            RowPart = 1;
         // printf("lineClusters: 3)hPad=%f RowPart=%f\n",hPad,RowPart);
         if (xloc1 < XPads) {
            iPad1 = 0;
         } else {
            iPad1 = Int_t((xloc1 - XPads) / fWpad);
         }
         if (xloc2 > -XPads) {
            iPad2 = nPads - 1;
         } else {
            iPad2 = Int_t((xloc2 - XPads) / fWpad);
         }
         // printf("lineClusters: irow=%d(%d-%d) XPads=%f FPart=%f LPart=%f \n",
         // irow,iPad1,iPad2,XPads,FirstPadPart,LastPadPart);
         for (Int_t ipad = iPad1; ipad <= iPad2; ipad++) {
            x1pad = XPads + ipad * fWpad;
            x2pad = x1pad + fWpad;
            if (xloc2 < x2pad) {
               FirstPadPart = 1.;
            } else {
               FirstPadPart = (x2pad - xloc1) / fWpad;
            }
            if (xloc1 > x1pad) {
               LastPadPart = 1.;
            } else {
               LastPadPart = (xloc2 - x1pad) / fWpad;
            }
            if (ipad == iPad1 && irow == iPad2) {
               PadPart = 1.;
            } else {
               if (irow == iRow1) {
                  PadPart = FirstPadPart;
               } else {
                  if (irow == iRow2) {
                     PadPart = LastPadPart;
                  } else {
                     PadPart = 1;
                  }
               }
            }
            // if(isec=1||isec==5)
            // printf("lineClusters: ipad=%d RowPart=%f PadPart=%f \n",ipad,RowPart,PadPart);
            SPads[ipad] = d_track * hPad * RowPart * PadPart;
         }
         nclust += ClustersInRow(irow, iPad1, iPad2, iRoc, SPads);
      }
      return nclust;
   }

   // - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - -

   Double_t s1, s2, b1, b2, x11, x12, x21, x22, tng;
   tng = fEsec.X() / fEsec.Y();
   for (Int_t irow = iRow1; irow <= iRow2; irow++) {
      b1         = fTpcSecGeo->LowRowEdge(irow);
      hPad       = fTpcSecGeo->PadHeight(irow);
      b2         = b1 + hPad;
      Int_t iRoc = fTpcSecGeo->NofRoc(irow);
      if (irow == iRow1 && irow == iRow2)
         RowPart = LowRowPart + UpRowPart - 1;
      else if (irow == iRow1)
         RowPart = LowRowPart;
      else if (irow == iRow2)
         RowPart = UpRowPart;
      else
         RowPart = 1;

      //  find min&max positions for the current row
      x11 = fR1sec.X() + (b1 - fR1sec.Y()) * tng;
      x12 = fR1sec.X() + (b2 - fR1sec.Y()) * tng;
      if (x11 < x12) {
         xloc1 = x11;
         xloc2 = x12;
      } else {
         xloc1 = x12;
         xloc2 = x11;
      }
      x21 = fR2sec.X() + (b1 - fR2sec.Y()) * tng;
      x22 = fR2sec.X() + (b2 - fR2sec.Y()) * tng;
      if (xloc1 > x21) xloc1 = x21;
      if (xloc2 < x21) xloc2 = x21;
      if (xloc1 > x22) xloc1 = x22;
      if (xloc2 < x22) xloc2 = x22;
      // if(xloc1>xloc2) printf("lineClusters: error\n");
      nPads = fTpcSecGeo->NofPadsInRow(irow);
      XPads = -nPads * fWpad / 2;
      if (xloc2 < XPads || xloc1 > -XPads) {
         continue;
      }
      // if(Nreconstructed==11&&isec==2)printf("lineClusters: x11=%f-(%f-%f)*%f=%f\n",fR1sec.X(),b1,fR1sec.Y(),tng,x11);
      // if(Nreconstructed==11&&isec==2)printf("lineClusters: x12=%f-(%f-%f)*%f=%f\n",fR1sec.X(),b2,fR1sec.Y(),tng,x12);
      // if(Nreconstructed==11&&isec==2)printf("lineClusters: x21=%f-(%f-%f)*%f=%f\n",fR2sec.X(),b1,fR2sec.Y(),tng,x21);
      // if(Nreconstructed==11&&isec==2)printf("lineClusters: x22=%f-(%f-%f)*%f=%f\n",fR2sec.X(),b2,fR2sec.Y(),tng,x22);
      if (xloc1 < XPads) {
         iPad1 = 0;
      } else {
         iPad1 = Int_t((xloc1 - XPads) / fWpad);
      }
      if (xloc2 > -XPads) {
         iPad2 = nPads - 1;
      } else { //
         if (xloc2 > 0) iPad2++;
         iPad2 = Int_t((xloc2 - XPads) / fWpad);
      }
      // if(Nreconstructed==11&&isec==2)printf("lineClusters: irow=%d nPads=%d XPads=%f xloc1,2(%f,%f)
      // iPad1,2(%d->%d)\n",irow,nPads,XPads,xloc1,xloc2,iPad1,iPad2);
      for (Int_t ipad = iPad1; ipad <= iPad2; ipad++) {
         s1          = XPads + ipad * fWpad;
         s2          = s1 + fWpad;
         spad        = padANDband(s1, s2, b1, b2);
         SPads[ipad] = spad * RowPart;
         // printf("lineClusters: %d-%d-%d S(%f,%f)=%f\n",isec,irow,ipad,(s1+s2)/2,(b1+b2)/2,SPads[ipad]);
         // if(Nreconstructed==11&&isec==2)printf("lineClusters: %d-%d-%d
         // (x11,x12,x21,x22)=(%f,%f,%f,%f)\n",isec,irow,ipad,x11,x12,x21,x22); printf("lineClusters: %d-%d-%d   S=%f
         // rP=%f)\n",isec,irow,ipad,spad,RowPart);
      }
      nclust += ClustersInRow(irow, iPad1, iPad2, iRoc, SPads);
   }
   return nclust;
}

//______________________________________________________________________

// find clusters in the row
Int_t TpcMiniMC::ClustersInRow(Int_t iRow, Int_t iPad1, Int_t iPad2, Int_t iRoc, Double_t *SPad)
{
   // if(Nreconstructed==10569)
   // printf("  ClustersInRow: s=%d %d(%d,%d)  S(%f ... %f)\n",
   // iSector,iRow,iPad1,iPad2,SPad[iPad1],SPad[iPad2]);
   //
   Double_t Signal, hpad;
   Int_t    nPadsInCluster = 0, nClust = 0;
   hpad = fTpcSecGeo->PadHeight(iRow);
   for (Int_t ipad = iPad1; ipad <= iPad2; ipad++) {
      // if(Nreconstructed==11&&iSector==2)printf("ClustersInRow: ipad=%d SPad=%f\n",ipad,SPad[ipad]);
      Signal = RandomPadSignal(iRoc, SPad[ipad]);
      if (Signal == 0) {
         if (nPadsInCluster > 0) {
            if (newCluster(hpad, nPadsInCluster)) nClust++; // if new cluster created
            nPadsInCluster = 0;
         }
         continue;
      } else {
         sPad[nPadsInCluster] = Signal;
         TVector2 LCxy        = fTpcSecGeo->LocalPadPosition(iRow, ipad);
         xPad[nPadsInCluster] = LCxy.X();
         yPad[nPadsInCluster] = LCxy.Y();
         if (!MFieldMode)
            zPad[nPadsInCluster] = -simLineZpos(LCxy);
         else
            zPad[nPadsInCluster] = -simHelixZpos(LCxy);
         nPadsInCluster++;
         if (nPadsInCluster == maxPcluster) {
            if (newCluster(hpad, nPadsInCluster)) nClust++; // if new cluster created
            // printf("ClustersInRow: ipad=%d nPadsInCluster=%d nClust=%d hits=%d\n",
            // ipad,nPadsInCluster,nClust,MpdTpcMiniMC_Nhits);
            nPadsInCluster = 0;
         }

         // if(Nreconstructed==11&&iSector==2)printf("ClustersInRow: event=%d iRow=%d ipad=%d signal=%f z=%f\n",
         // Nreconstructed,iRow,ipad,Signal,zPad[nPadsInCluster]);
      }
   }
   if (nPadsInCluster > 0) {
      if (newCluster(hpad, nPadsInCluster)) nClust++; // if new cluster created
   }
   return nClust;
}

//______________________________________________________________________

// create new cluster
Bool_t TpcMiniMC::newCluster(Double_t hpad, Int_t nPads)
{
   Double_t SPAD = 0., xClaster = 0., yClaster = 0., zClaster = 0.;
   Double_t cosAB, Zalong, sigma;
   TVector3 eef, A, B, AmB, amb, rlClaster, rgClaster, rbClaster;
   for (Int_t ipad = 0; ipad < nPads; ipad++) {
      SPAD += sPad[ipad];
      xClaster += xPad[ipad] * sPad[ipad];
      yClaster += yPad[ipad] * sPad[ipad];
      zClaster += zPad[ipad] * sPad[ipad];
      // printf("newCluster: ipad=%d rL(%f,%f,%f)\n",ipad,xPad[ipad],yPad[ipad],zPad[ipad]);
   }
   xClaster /= SPAD;
   yClaster /= SPAD;
   zClaster /= SPAD;

   /*
      <rgClaster> - global coordinates of the hit
      <rlClaster> - local coordinates of the hit
      <e_3> - unit vector in GSC along the Z-axis of LSC
      <eef> - unit vector in GSC of the electric field
      cosAB=(e_3*eef) <A>=Zalong*eef  <B>=(-Zalong*cosAB)*<e_3>
      <AmB>=<A>-<Bn>  <amb>={glo2loc}*<AmB>
         xl=Xl-<amb>_x  yl=Yl-<amb>_y  zln=--Zalong*cosAB
         <Xg,Yg,Zg>=R0shift+{loc2glo}*<xln,Yln,zln=-Zalong*cosAB>
   */
   if (e_3.Z() > 0)
      eef = eEf.Unit();
   else
      eef = -eEf.Unit();
   cosAB     = e_3 * eef;
   Zalong    = TMath::Abs(zClaster);
   A         = Zalong * eef;
   B         = (Zalong * cosAB) * e_3;
   AmB       = A - B;
   amb       = glo2loc * AmB;
   rlClaster = TVector3(xClaster - amb.X(), yClaster - amb.Y(), -Zalong * cosAB);
   rgClaster = R0shift + loc2glo * rlClaster;
   if (rgClaster.Z() < fZmin || rgClaster.Z() > fZmax) return false;
   // printf("newCluster: rG(%f,%f,%f)\n",rgClaster.X(),rgClaster.Y(),rgClaster.Z());
   sigma = TMath::Sqrt((fWpad * fWpad + hpad * hpad + (Zalong * sig_Z_pos) * (Zalong * sig_Z_pos)) / nPads);
   //    sigma=1./SPAD;
   MpdTpcMiniMC_WHits[MpdTpcMiniMC_Nhits] = sigma;
   MpdTpcMiniMC_iSec[MpdTpcMiniMC_Nhits]  = iSector;
   MpdTpcMiniMC_XHits[MpdTpcMiniMC_Nhits] = rgClaster.X();
   MpdTpcMiniMC_YHits[MpdTpcMiniMC_Nhits] = rgClaster.Y();
   MpdTpcMiniMC_ZHits[MpdTpcMiniMC_Nhits] = rgClaster.Z();
   if (MFieldMode) {
      rbClaster                               = G2B * rgClaster;
      MpdTpcMiniMC_XBHits[MpdTpcMiniMC_Nhits] = rbClaster.X();
      MpdTpcMiniMC_YBHits[MpdTpcMiniMC_Nhits] = rbClaster.Y();
      MpdTpcMiniMC_ZBHits[MpdTpcMiniMC_Nhits] = rbClaster.Z();
   }
   MpdTpcMiniMC_Nhits++;
   // printf("%dth newCluster of %dpads rL(%f,%f,%f) rG(%f,%f,%f)\n",
   // MpdTpcMiniMC_Nhits,nPads,rlClaster.X(),rlClaster.Y(),rlClaster.Z(),
   // rgClaster.X(),rgClaster.Y(),rgClaster.Z());

   return true;
}

//______________________________________________________________________

// simulate z pozition
Double_t TpcMiniMC::simLineZpos(TVector2 LCxy)
{

   // a) lower the perpendicular from LCxy to the track center line
   // b) find the distance from the base of the perpendicular along the vector E
   //    to the track of the particle.
   // This will be the average value of the Z coordinate measurement
   Double_t t    = ((LCxy.X() - fR0sec.X()) * uEsec.X() + (LCxy.Y() - fR0sec.Y()) * uEsec.Y());
   TVector3 SecL = TVector3(LCxy.X(), LCxy.Y(), 0.);
   TVector3 SecG = R0shift + loc2glo * SecL;
   //  TVector3 baseL=TVector3(fR0sec.X()+fEsec.X()*t,fR0sec.Y()+fEsec.Y()*t,0.);
   TVector3 baseL = TVector3(fR0sec.X() + uEsec.X() * t, fR0sec.Y() + uEsec.Y() * t, 0.);
   TVector3 baseG = R0shift + loc2glo * baseL;
   if (TMath::Abs(fEsec.X()) > epsMax) {
      t = (baseL.X() - fR0sec.X()) / fEsec.X();
   } else {
      t = (baseL.Y() - fR0sec.Y()) / fEsec.Y();
   }
   TVector3 base2track = fR1_c + (fR2_c - fR1_c) * t; // corresponding track point
   TVector3 D          = base2track - SecG;           // distance to the track point from the cluster
   D                   = base2track - baseG;          // distance to the track point from the cluster
   Double_t z          = fRandom->Gaus(D.Mag(), D.Mag() * sig_Z_pos);
   return z;
}

//______________________________________________________________________

// find particle cross points (no MF)
void TpcMiniMC::crossLinePoints()
{
   // 1) the track can cross chamber bounds and the membrana
   //    in this case points fP11,fP12 corresponds to the
   //    chamber #1 (z>0) and points fP21,fP22 to the chamber #2 (z<0)
   // 2) The track does not crosses the sensitive volume of a chamber,
   //    in this case |fP11| is zero for the chamber #1 or |fP21|=0
   //    for the chamber #2

   Double_t a, b, c, t_0, t_1, t_2, z_1, z_2, az1, az2;
   // find particle cross points with the chamber cilinder
   a = pLas.X() * pLas.X() + pLas.Y() * pLas.Y();
   if (TMath::Sqrt(a) > epsMin) { // pLas is no parallel to Z-axis
      b   = pLas.X() * r0Beam.X() + pLas.Y() * r0Beam.Y();
      c   = r0Beam.X() * r0Beam.X() + r0Beam.Y() * r0Beam.Y() - fr_svol * fr_svol;
      t_1 = (-b + sqrt(b * b - a * c)) / a;
      z_1 = r0Beam.Z() + pLas.Z() * t_1;
      t_2 = (-b - sqrt(b * b - a * c)) / a;
      z_2 = r0Beam.Z() + pLas.Z() * t_2;
      if (z_1 > z_2) {
         Double_t w = z_2;
         z_2        = z_1;
         z_1        = w;
         w          = t_2;
         t_2        = t_1;
         t_1        = w;
      }
      if (z_1 < -fZmax) { // z_1 is behind the 2nd chamber PP
         t_1 = (-fZmax - r0Beam.Z()) / pLas.Z();
         z_1 = -fZmax;
      }
      if (z_2 > fZmax) { // z_1 is behind the 2nd chamber PP
         t_2 = (fZmax - r0Beam.Z()) / pLas.Z();
         z_2 = fZmax;
      }
      // printf("crossLinePoints: t1,2(%f,%f) z_1,2(%f,%f)\n",t_1,t_2,z_1,z_2);
      if (z_1 * z_2 < 0) { // track crosses both chambers
         if (z_1 > -fZmin) {
            fP21.SetXYZ(0., 0., 0.);
            if (z_2 < fZmin)
               fP11.SetXYZ(0., 0., 0.);
            else {
               fP11 = r0Beam + pLas * (fZmin - r0Beam.Z() / pLas.Z());
               fP12 = r0Beam + pLas * t_2;
            }
         } else {
            fP21 = r0Beam + pLas * t_1;
            fP22 = r0Beam + pLas * ((-fZmin - r0Beam.Z()) / pLas.Z());
            if (z_2 < fZmin) {
               fP11.SetXYZ(0., 0., 0.);
            } else {
               fP11 = r0Beam + pLas * ((fZmin - r0Beam.Z()) / pLas.Z());
               fP12 = r0Beam + pLas * t_2;
            }
         }
      } else {
         if (z_1 >= 0) {
            fP21.SetXYZ(0., 0., 0.);
            if (z_1 < fZmin) {
               if (z_2 > fZmin) {
                  fP11 = r0Beam + pLas * ((fZmin - r0Beam.Z()) / pLas.Z());
                  fP12 = r0Beam + pLas * t_2;
               } else
                  fP11.SetXYZ(0., 0., 0.);
            } else {
               fP11 = r0Beam + pLas * t_1;
               fP12 = r0Beam + pLas * t_2;
            }
         } else {
            fP11.SetXYZ(0., 0., 0.);
            if (z_2 > -fZmin) {
               if (z_1 < -fZmin) {
                  fP22 = r0Beam + pLas * ((-fZmin - r0Beam.Z()) / pLas.Z());
                  fP12 = r0Beam + pLas * t_1;
               } else
                  fP21.SetXYZ(0., 0., 0.);
            } else {
               fP21 = r0Beam + pLas * t_1;
               fP22 = r0Beam + pLas * t_2;
            }
         }
      }
   } else {
      fP11.SetXYZ(0., 0., 0.);
      fP21.SetXYZ(0., 0., 0.);
   }
   // printf("crossLinePoints: pLas_in(%f,%f,%f)\n",pLas.X(),pLas.Y(),pLas.Z());
   // printf("crossLinePoints: fP11(%f,%f,%f)   fP12(%f,%f,%f)\n",
   // fP11.X(),fP11.Y(),fP11.Z(),fP12.X(),fP12.Y(),fP12.Z());
   // printf("crossLinePoints: fP21(%f,%f,%f)   fP22(%f,%f,%f)\n",
   // fP21.X(),fP21.Y(),fP21.Z(),fP22.X(),fP22.Y(),fP22.Z());
}

//______________________________________________________________________

// find the strait band on the sector plane, false if doesn't cross
Bool_t TpcMiniMC::crossCheck(Int_t iSec, TVector3 Tr1g, TVector3 Tr2g)
{
   // Tr1g - point 1 on the track (GlobalCS)
   // Tr2g - point 2 on the track (GlobalCS)
   // the function is "true" if the projection of the track interval
   // [Tr1g,Tr2g] on the sector #iSec (interval [xg_l1,xg_l2] in LSC)
   // crosses any pad of the sector.
   // In this case in LCS
   // fR0sec - the point on the track line
   // fR1sec,fR2sec - points on the left and right track boundaries
   // fEsec - the vector from xg_l1 xg_l2 projections of Tr1g,Tr2g jn pads plane

   // if(Nreconstructed==10569) printf("crossCheck: sec=%d T1(%f,%f,%f) T2(%f,%f,%f)\n",
   // iSec,Tr1g.X(),Tr1g.Y(),Tr1g.X(),Tr2g.X(),Tr2g.Y(),Tr2g.Z());
   const Int_t isec = iSec;
   // Get loc2glo transformation parameters for the sector isec
   fTpcSecGeo->SectorTransformation(isec, R0shift, loc2glo);
   e_3            = TVector3(loc2glo.XZ(), loc2glo.YZ(), loc2glo.ZZ());
   Double_t t_1   = (e_3 * (R0shift - Tr1g)) / (e_3 * eEf); // to 1st sector plane point
   Double_t t_2   = (e_3 * (R0shift - Tr2g)) / (e_3 * eEf); // to 2nd sector plane point
   TVector3 xg_l1 = Tr1g + t_1 * eEf;                       // Point 1 on the sector plane (GCS)
   TVector3 xg_l2 = Tr2g + t_2 * eEf;                       // Point 2 on the sector plane (GCS)
   fTpcSecGeo->SectorBackTransformation(isec, G0shift, glo2loc);
   TVector3 xl_l1 = G0shift + glo2loc * xg_l1;      // Point 1 on the sector plane (LCS)
   fR0sec         = TVector2(xl_l1.X(), xl_l1.Y()); // Point on the track center line
   TVector3 xl_l2 = G0shift + glo2loc * xg_l2;      // Point 2 on the track center line
   fR0sec1        = TVector2(xl_l2.X(), xl_l2.Y());
   fEsec          = fR0sec1 - fR0sec;
   uEsec          = fEsec / fEsec.Mod();
   // Check - does it cross sector pads?
   Double_t y_max = fTpcSecGeo->GetMaxYl();
   Double_t y_min = fTpcSecGeo->GetMinYl();
   Double_t t, x11, x12, x21, x22, xy_max, xy_min, y1min, y2min, y1max, y2max;
   t          = 1. / fEsec.Mod();
   Bool_t Yes = false; // true, if a particle leaves a track on pads area
   if (TMath::Abs(fEsec.Y()) > epsMin) {
      if (fEsec.Y() < 0)
         Uniper = TVector2(t * fEsec.Y(), -t * fEsec.X());
      else
         Uniper = TVector2(-t * fEsec.Y(), t * fEsec.X());
      fR1sec  = fR0sec + dh_track * Uniper;
      fR1sec1 = fR1sec + fEsec;
      fR2sec  = fR0sec - dh_track * Uniper;
      fR2sec1 = fR2sec + fEsec;
      if (fR1sec.Y() < fR1sec1.Y()) {
         y1min = fR1sec.Y();
         y1max = fR1sec1.Y();
      } else {
         y1min = fR1sec1.Y();
         y1max = fR1sec.Y();
      }
      if (fR2sec.Y() < fR2sec1.Y()) {
         y2min = fR2sec.Y();
         y2max = fR2sec1.Y();
      } else {
         y2min = fR2sec1.Y();
         y2max = fR2sec.Y();
      }
      if ((y1max <= y_min && y2max <= y_min) || (y1min >= y_max && y2min >= y_min)) return Yes;
      // a) top bound
      t      = (y_max - fR1sec.Y()) / fEsec.Y();
      x12    = fR1sec.X() + fEsec.X() * t;
      t      = (y_max - fR2sec.Y()) / fEsec.Y();
      x22    = fR2sec.X() + fEsec.X() * t;
      xy_max = fTpcSecGeo->GetMaxX();
      if (x12 >= x22) printf("CrossCheck: algorithm error1!\n");
      if ((-xy_max < x12 && x12 < xy_max) || (-xy_max < x22 && x22 < xy_max)) {
         Yes = true;
      } else {
         // b) bottom bound
         t   = (y_min - fR1sec.Y()) / fEsec.Y();
         x11 = fR1sec.X() + fEsec.X() * t;
         t   = (y_min - fR2sec.Y()) / fEsec.Y();
         x21 = fR2sec.X() + fEsec.X() * t;
         if (x11 >= x21) printf("CrossCheck: algorithm error2!\n");
         if ((-xy_min < x11 && x11 < xy_min) || (-xy_min < x21 && x21 < xy_min)) {
            Yes = true;
         } else {
            // c) check the simmetry line cross
            if ((x22 <= -xy_max && x11 >= xy_min) || (x12 >= xy_max && x21 <= -xy_min)) {
               Yes = true;
            }
         }
      }
   } else {
      if (y_min - dh_track <= fR0sec.Y() && fR0sec.Y() <= y_max + dh_track) {
         Yes     = true;
         fR1sec  = TVector2(fR0sec.X(), fR0sec.Y() + dh_track);
         fR1sec1 = fR1sec + fEsec;
         fR2sec  = TVector2(fR0sec.X(), fR0sec.Y() - dh_track);
         fR2sec1 = fR2sec + fEsec;
      }
   }
   //    if(Yes) { // find the band from the direct line track
   fR1_c = Tr1g;
   fR2_c = Tr2g;
   //    }
   // printf("crossCheck: sector=%d  Yes=%d\n",iSec,Yes);
   return Yes;
}

//______________________________________________________________________
// simulate theta angle for cosmic muon
Double_t TpcMiniMC::theta_cosm()
{
   Double_t q = fRandom->Rndm();
   if (k_par == 1) {
      return 0.5 * TMath::Pi() * (sqrt(2) - sqrt(2 - 2 * q));
   } else {
      if (k_par == 2) {
         return TMath::ACos(q);
      } else {
         cout << " MiniMC: k_par is out of range: " << k_par << endl;
         return 0;
      }
   }
}

//______________________________________________________________________
// simulate theta angle for beam
Double_t TpcMiniMC::theta_beam()
{
   Double_t q = fRandom->Rndm();
   return TMath::ASin(q);
}

//______________________________________________________________________
// simulate pad signal
Double_t TpcMiniMC::RandomPadSignal(Int_t iRoc, Double_t SPad)
{
   Double_t sigma     = SPad * sig_padS[iRoc];
   Double_t PadSignal = fRandom->Gaus(SPad, sigma);
   if (PadSignal < 0.) {
      return 0.;
   } else {
      if (PadSignal > fTpcSecGeo->GetPadSforROC(iRoc)) {
         return fTpcSecGeo->GetPadSforROC(iRoc);
      } else {
         return PadSignal;
      }
   }
}

//______________________________________________________________________
// calculate the pad square covered by the track band
Double_t TpcMiniMC::padANDband(Double_t s1, Double_t s2, Double_t b1, Double_t b2)
{

   Int_t    nP1 = 0, nM1 = 0, nP2 = 0, nM2 = 0, error;
   Bool_t   Vl1pos[4], Vl2pos[4];
   Double_t spad, SPad, xp1, yp1, xp2, yp2, stra, stri;
   SPad  = (b2 - b1) * (s2 - s1);
   error = 0;
   // calculate number of positive & negative vertexes for band bounderies
   if (f1(s1, b1) > 0) {
      nP1++;
      Vl1pos[0] = true;
   } else {
      nM1++;
      Vl1pos[0] = false;
   }
   if (f1(s1, b2) > 0) {
      nP1++;
      Vl1pos[1] = true;
   } else {
      nM1++;
      Vl1pos[1] = false;
   }
   if (f1(s2, b2) > 0) {
      nP1++;
      Vl1pos[2] = true;
   } else {
      nM1++;
      Vl1pos[2] = false;
   }
   if (f1(s2, b1) > 0) {
      nP1++;
      Vl1pos[3] = true;
   } else {
      nM1++;
      Vl1pos[3] = false;
   }

   if (f2(s1, b1) > 0) {
      nP2++;
      Vl2pos[0] = true;
   } else {
      nM2++;
      Vl2pos[0] = false;
   }
   if (f2(s1, b2) > 0) {
      nP2++;
      Vl2pos[1] = true;
   } else {
      nM2++;
      Vl2pos[1] = false;
   }
   if (f2(s2, b2) > 0) {
      nP2++;
      Vl2pos[2] = true;
   } else {
      nM2++;
      Vl2pos[2] = false;
   }
   if (f2(s2, b1) > 0) {
      nP2++;
      Vl2pos[3] = true;
   } else {
      nM2++;
      Vl2pos[3] = false;
   }
   // if(Nreconstructed==10569) printf("padANDband: Vl1pos=(%d,%d,%d,%d) Vl2pos=(%d,%d,%d,%d) nP1=%d nM1=%d nP2=%d
   // nM2=%d\n", Vl1pos[0],Vl1pos[1],Vl1pos[2],Vl1pos[3],Vl2pos[0],Vl2pos[1],Vl2pos[2],Vl2pos[3],nP1,nM1,nP2,nM2);
   //  a) full pad
   if ((nM1 == 4 && nP2 == 4) || (nP1 == 4 && nM2 == 4)) return SPad;
   // b) pad minus a triangle
   if (nM1 == 1 && nM2 == 4) {
      if (!Vl1pos[0]) {
         return SPad - sTriangle(1, s1, b1);
      } else {
         if (!Vl1pos[1]) {
            return SPad - sTriangle(1, s1, b2);
         } else {
            if (!Vl1pos[2]) {
               return SPad - sTriangle(1, s2, b2);
            } else {
               if (!Vl1pos[3]) {
                  return SPad - sTriangle(1, s2, b1);
               } else {
                  cout << "MpdTpcMiniMC: padANDband's error1" << endl;
                  return -5;
               }
            }
         }
      }
   }

   if (nP1 == 1 && nP2 == 4) {
      if (Vl1pos[0]) {
         return SPad - sTriangle(1, s1, b1);
      } else {
         if (Vl1pos[1]) {
            return SPad - sTriangle(1, s1, b2);
         } else {
            if (Vl1pos[2]) {
               return SPad - sTriangle(1, s2, b2);
            } else {
               if (Vl1pos[3]) {
                  return SPad - sTriangle(1, s2, b1);
               } else {
                  cout << "MpdTpcMiniMC: padANDband's error2" << endl;
                  return -5;
               }
            }
         }
      }
   }

   if (nM2 == 1 && nM1 == 4) {
      if (!Vl2pos[0]) {
         return SPad - sTriangle(2, s1, b1);
      } else {
         if (!Vl2pos[1]) {
            return SPad - sTriangle(2, s1, b2);
         } else {
            if (!Vl2pos[2]) {
               return SPad - sTriangle(2, s2, b2);
            } else {
               if (!Vl2pos[3]) {
                  return SPad - sTriangle(2, s2, b1);
               } else {
                  cout << "MpdTpcMiniMC: padANDband's error3" << endl;
                  return -5;
               }
            }
         }
      }
   }

   if (nP2 == 1 && nP1 == 4) {
      if (Vl2pos[0]) {
         return SPad - sTriangle(2, s1, b1);
      } else {
         if (Vl2pos[1]) {
            return SPad - sTriangle(2, s1, b2);
         } else {
            if (Vl2pos[2]) {
               return SPad - sTriangle(2, s2, b2);
            } else {
               if (Vl2pos[3]) {
                  return SPad - sTriangle(2, s2, b1);
               } else {
                  cout << "MpdTpcMiniMC: padANDband's error4" << endl;
                  return -5;
               }
            }
         }
      }
   }

   // c) pad minus a trapezoid
   if (nM1 == 2 && nP2 == 4) {
      if (!(Vl1pos[0] || Vl1pos[1])) {
         return sTrapezoid(1, 1, s1, b1, b2);
      } else {
         if (!(Vl1pos[1] || Vl1pos[2])) {
            return sTrapezoid(1, 2, b2, s1, s2);
         } else {
            if (!(Vl1pos[2] || Vl1pos[3])) {
               return sTrapezoid(1, 1, s2, b2, b1);
            } else {
               if (!(Vl1pos[3] || Vl1pos[0])) {
                  return sTrapezoid(1, 2, b1, s2, s1);
               } else {
                  cout << "MpdTpcMiniMC: padANDband's error5" << endl;
                  return -5;
               }
            }
         }
      }
   }

   if (nP1 == 2 && nM2 == 4) {
      if (Vl1pos[0] && Vl1pos[1]) {
         return sTrapezoid(1, 1, s1, b1, b2);
      } else {
         if (Vl1pos[1] && Vl1pos[2]) {
            return sTrapezoid(1, 2, b2, s1, s2);
         } else {
            if (Vl1pos[2] && Vl1pos[3]) {
               return sTrapezoid(1, 1, s2, b2, b1);
            } else {
               if (Vl1pos[3] && Vl1pos[0]) {
                  return sTrapezoid(1, 2, b1, s2, s1);
               } else {
                  cout << "MpdTpcMiniMC: padANDband's error6" << endl;
                  return -5;
               }
            }
         }
      }
   }

   if (nM2 == 2 && nP1 == 4) {
      if (!(Vl2pos[0] || Vl2pos[1])) {
         return sTrapezoid(2, 1, s1, b1, b2);
      } else {
         if (!(Vl2pos[1] || Vl2pos[2])) {
            return sTrapezoid(2, 2, b2, s1, s2);
         } else {
            if (!(Vl2pos[2] || Vl2pos[3])) {
               return sTrapezoid(2, 1, s2, b2, b1);
            } else {
               if (!(Vl2pos[3] || Vl2pos[0])) {
                  return sTrapezoid(2, 2, b1, s2, s1);
               } else {
                  cout << "MpdTpcMiniMC: padANDband's error7" << endl;
                  return -5;
               }
            }
         }
      }
   }

   if (nP2 == 2 && nM1 == 4) {
      if (Vl2pos[0] && Vl2pos[1]) {
         return sTrapezoid(2, 1, s1, b1, b2);
      } else {
         if (Vl2pos[1] && Vl2pos[2]) {
            return sTrapezoid(2, 2, b2, s1, s2);
         } else {
            if (Vl2pos[2] && Vl2pos[3]) {
               return sTrapezoid(2, 1, s2, b2, b1);
            } else {
               if (Vl2pos[3] && Vl2pos[0]) {
                  return sTrapezoid(2, 2, b1, s2, s1);
               } else {
                  cout << "MpdTpcMiniMC: padANDband's error8" << endl;
                  return -5;
               }
            }
         }
      }
   }

   // d) triangle
   if (nM1 == 1 && nP2 == 4) {
      if (!Vl1pos[0]) {
         return sTriangle(1, s1, b1);
      } else {
         if (!Vl1pos[1]) {
            return sTriangle(1, s1, b2);
         } else {
            if (!Vl1pos[2]) {
               return sTriangle(1, s2, b2);
            } else {
               if (!Vl1pos[3]) {
                  return sTriangle(1, s2, b1);
               } else {
                  cout << "MpdTpcMiniMC: padANDband's error9" << endl;
                  return -5;
               }
            }
         }
      }
   }

   if (nP1 == 1 && nM2 == 4) {
      if (Vl1pos[0]) {
         return sTriangle(1, s1, b1);
      } else {
         if (Vl1pos[1]) {
            return sTriangle(1, s1, b2);
         } else {
            if (Vl1pos[2]) {
               return sTriangle(1, s2, b2);
            } else {
               if (Vl1pos[3]) {
                  return sTriangle(1, s2, b1);
               } else {
                  cout << "MpdTpcMiniMC: padANDband's error10" << endl;
                  return -5;
               }
            }
         }
      }
   }

   if (nM2 == 1 && nP1 == 4) {
      if (!Vl2pos[0]) {
         return sTriangle(2, s1, b1);
      } else {
         if (!Vl2pos[1]) {
            return sTriangle(2, s1, b2);
         } else {
            if (!Vl2pos[2]) {
               return sTriangle(2, s2, b2);
            } else {
               if (!Vl2pos[3]) {
                  return sTriangle(2, s2, b1);
               } else {
                  cout << "MpdTpcMiniMC: padANDband's error11" << endl;
                  return -5;
               }
            }
         }
      }
   }

   if (nP2 == 1 && nM1 == 4) {
      if (Vl2pos[0]) {
         return sTriangle(2, s1, b1);
      } else {
         if (Vl2pos[1]) {
            return sTriangle(2, s1, b2);
         } else {
            if (Vl2pos[2]) {
               return sTriangle(2, s2, b2);
            } else {
               if (Vl2pos[3]) {
                  return sTriangle(2, s2, b1);
               } else {
                  cout << "MpdTpcMiniMC: padANDband's error12" << endl;
                  return -5;
               }
            }
         }
      }
   }

   // e) pad minus 2 triangles
   if (nM1 == 1 && nP2 == 1) {
      if (!Vl1pos[0] && Vl2pos[2]) {
         return SPad - sTriangle(1, s1, b1) - sTriangle(2, s2, b2);
      } else {
         if (!Vl1pos[1] && Vl2pos[3]) {
            return SPad - sTriangle(1, s1, b2) - sTriangle(2, s2, b1);
         } else {
            if (!Vl1pos[2] && Vl2pos[0]) {
               return SPad - sTriangle(1, s2, b2) - sTriangle(2, s1, b1);
            } else {
               if (!Vl1pos[3] && Vl2pos[1]) {
                  return SPad - sTriangle(1, s2, b1) - sTriangle(2, s1, b2);
               } else {
                  cout << "MpdTpcMiniMC: padANDband's error13" << endl;
                  return -5;
               }
            }
         }
      }
   }

   if (nP1 == 1 && nM2 == 1) {
      if (Vl1pos[0] && !Vl2pos[2]) {
         return SPad - sTriangle(1, s1, b1) - sTriangle(2, s2, b2);
      } else {
         if (Vl1pos[1] && !Vl2pos[3]) {
            return SPad - sTriangle(1, s1, b2) - sTriangle(2, s2, b1);
         } else {
            if (Vl1pos[2] && !Vl2pos[0]) {
               return SPad - sTriangle(1, s2, b2) - sTriangle(2, s1, b1);
            } else {
               if (Vl1pos[3] && !Vl2pos[1]) {
                  return SPad - sTriangle(1, s2, b1) - sTriangle(2, s1, b2);
               } else {
                  cout << "MpdTpcMiniMC: padANDband's error14" << endl;
                  return -5;
               }
            }
         }
      }
   }

   // f) trapezoid minus triangle (case 6)
   if (nM1 == 2 && nM2 == 1) {
      if (!Vl2pos[0]) {
         if (!Vl1pos[0] && !Vl1pos[3])
            return sTrapezoid(1, 2, b1, s1, s2) - sTriangle(2, s1, b1);
         else if (!Vl1pos[0] && !Vl1pos[1])
            return sTrapezoid(1, 1, s1, b1, b2) - sTriangle(2, s1, b1);
      } else if (!Vl2pos[1]) {
         if (!Vl1pos[1] && !Vl1pos[2])
            return sTrapezoid(1, 2, b2, s1, s2) - sTriangle(2, s1, b2);
         else if (!Vl1pos[1] && !Vl1pos[0])
            return sTrapezoid(1, 1, s1, b1, b2) - sTriangle(2, s1, b2);
      } else if (!Vl2pos[2]) {
         if (!Vl1pos[2] && !Vl1pos[3])
            return sTrapezoid(1, 1, s2, b1, b2) - sTriangle(2, s2, b2);
         else if (!Vl1pos[2] && !Vl1pos[1]) {
            return sTrapezoid(1, 2, b2, s1, s2) - sTriangle(2, s2, b2);
         }
      } else if (!Vl2pos[3]) {
         if (!Vl1pos[3] && !Vl1pos[0])
            return sTrapezoid(1, 2, b1, s1, s2) - sTriangle(2, s2, b1);
         else if (!Vl1pos[3] && !Vl1pos[2])
            return sTrapezoid(1, 1, s2, b1, b2) - sTriangle(2, s2, b1);
      }
      cout << "MpdTpcMiniMC: padANDband's error61" << endl;
      exit(61);
   }
   // - - - - - - - -
   if (nM1 == 1 && nM2 == 2) {
      if (!Vl1pos[0]) {
         if (!Vl2pos[0] && !Vl2pos[3])
            return sTrapezoid(2, 2, b1, s1, s2) - sTriangle(1, s1, b1);
         else if (!Vl2pos[0] && !Vl2pos[1])
            return sTrapezoid(2, 1, s1, b1, b2) - sTriangle(1, s1, b1);
      } else if (!Vl1pos[1]) {
         if (!Vl2pos[1] && !Vl2pos[2])
            return sTrapezoid(2, 2, b2, s1, s2) - sTriangle(1, s1, b2);
         else if (!Vl2pos[1] && !Vl2pos[0])
            return sTrapezoid(2, 1, s1, b1, b2) - sTriangle(1, s1, b2);
      } else if (!Vl1pos[2]) {
         if (!Vl2pos[2] && !Vl2pos[3])
            return sTrapezoid(2, 1, s2, b1, b2) - sTriangle(1, s2, b2);
         else if (!Vl2pos[2] && !Vl2pos[1])
            return sTrapezoid(2, 2, b2, s1, s2) - sTriangle(1, s2, b2);
      } else if (!Vl1pos[3]) {
         if (!Vl2pos[3] && !Vl2pos[0])
            return sTrapezoid(2, 2, b1, s1, s2) - sTriangle(1, s2, b1);
         else if (!Vl2pos[3] && !Vl2pos[2])
            return sTrapezoid(2, 1, s2, b1, b2) - sTriangle(1, s2, b1);
      }
      cout << "MpdTpcMiniMC: padANDband's error62" << endl;
      exit(62);
   }
   // - - - - - - - -
   if (nP1 == 2 && nP2 == 1) {
      if (Vl2pos[0]) {
         if (Vl1pos[0] && Vl1pos[3])
            return sTrapezoid(1, 2, b1, s1, s2) - sTriangle(2, s1, b1);
         else if (Vl1pos[0] && Vl1pos[1])
            return sTrapezoid(1, 1, s1, b1, b2) - sTriangle(2, s1, b1);
      } else if (Vl2pos[1]) {
         if (Vl1pos[1] && Vl1pos[2])
            return sTrapezoid(1, 2, b2, s1, s2) - sTriangle(2, s1, b2);
         else if (Vl1pos[1] && Vl1pos[0])
            return sTrapezoid(1, 1, s1, b1, b2) - sTriangle(2, s1, b2);
      } else if (Vl2pos[2]) {
         if (Vl1pos[2] && Vl1pos[3])
            return sTrapezoid(1, 1, s2, b1, b2) - sTriangle(2, s2, b2);
         else if (Vl1pos[2] && Vl1pos[1])
            return sTrapezoid(1, 2, b2, s1, s2) - sTriangle(2, s2, b2);
      } else if (Vl2pos[3]) {
         if (Vl1pos[3] && Vl1pos[0])
            return sTrapezoid(1, 2, b1, s1, s2) - sTriangle(2, s2, b1);
         else if (Vl1pos[3] && Vl1pos[2])
            return sTrapezoid(1, 1, s2, b1, b2) - sTriangle(2, s2, b1);
      }
      cout << "MpdTpcMiniMC: padANDband's error63" << endl;
      exit(63);
   }
   // - - - - - - - -
   if (nP1 == 1 && nP2 == 2) {
      if (Vl1pos[0]) {
         if (Vl2pos[0] && Vl2pos[3])
            return sTrapezoid(2, 2, b1, s1, s2) - sTriangle(1, s1, b1);
         else if (Vl2pos[0] && Vl2pos[1])
            return sTrapezoid(2, 1, s1, b1, b2) - sTriangle(1, s1, b1);
      } else if (Vl1pos[1]) {
         if (Vl2pos[1] && Vl2pos[2])
            return sTrapezoid(2, 2, b2, s1, s2) - sTriangle(1, s1, b2);
         else if (Vl2pos[1] && Vl2pos[0])
            return sTrapezoid(2, 1, s1, b1, b2) - sTriangle(1, s1, b2);
      } else if (Vl1pos[2]) {
         if (Vl2pos[2] && Vl2pos[3])
            return sTrapezoid(2, 1, s2, b1, b2) - sTriangle(1, s2, b2);
         else if (Vl2pos[2] && Vl2pos[1])
            return sTrapezoid(2, 2, b2, s1, s2) - sTriangle(1, s2, b2);
      } else if (Vl1pos[3]) {
         if (Vl2pos[3] && Vl2pos[0])
            return sTrapezoid(2, 2, b1, s1, s2) - sTriangle(1, s2, b1);
         else if (Vl2pos[3] && Vl2pos[2])
            return sTrapezoid(2, 1, s2, b1, b2) - sTriangle(1, s2, b1);
      }
      cout << "MpdTpcMiniMC: padANDband's error64" << endl;
      exit(64);
   }

   // g) triangle minus triangle (case 7)
   if (nM1 == 1 && nM2 == 1) {
      if (!Vl1pos[0] && !Vl2pos[0])
         return TMath::Abs(sTriangle(1, s1, b1) - sTriangle(2, s1, b1));
      else if (!Vl1pos[1] && !Vl2pos[1])
         return TMath::Abs(sTriangle(1, s1, b2) - sTriangle(2, s1, b2));
      else if (!Vl1pos[2] && !Vl2pos[2])
         return TMath::Abs(sTriangle(1, s2, b2) - sTriangle(2, s2, b2));
      else if (!Vl1pos[3] && !Vl2pos[3])
         return TMath::Abs(sTriangle(1, s2, b1) - sTriangle(2, s2, b1));
      else {
         cout << "MpdTpcMiniMC: padANDband's error71" << endl;
         exit(71);
      }
   }
   // - - - - - - - -
   if (nP1 == 1 && nP2 == 1) {
      if (Vl1pos[0] && Vl2pos[0])
         return TMath::Abs(sTriangle(1, s1, b1) - sTriangle(2, s1, b1));
      else if (Vl1pos[1] && Vl2pos[1])
         return TMath::Abs(sTriangle(1, s1, b2) - sTriangle(2, s1, b2));
      else if (Vl1pos[2] && Vl2pos[2])
         return TMath::Abs(sTriangle(1, s2, b2) - sTriangle(2, s2, b2));
      else if (Vl1pos[3] && Vl2pos[3])
         return TMath::Abs(sTriangle(1, s2, b1) - sTriangle(2, s2, b1));
      else {
         cout << "MpdTpcMiniMC: padANDband's error72" << endl;
         exit(72);
      }
   }

   // h) parallelogram (case 8)

   if (nP1 == 2 && nP2 == 2) {
      if ((Vl1pos[0] && Vl1pos[1]) || (Vl1pos[2] && Vl1pos[3]))
         return sParallelogram(2, s1, s2, b1, b2);
      else if ((Vl1pos[1] && Vl1pos[2]) || (Vl1pos[3] && Vl1pos[0])) {
         // if(Nreconstructed==10569) printf("padANDband: s1,s2,b1,b2=(%f,%f,%f,%f)\n",s1,s2,b1,b2);
         // if(Nreconstructed==10569) printf("padANDband: %f\n",sParallelogram(2,s1,s2,b1,b2));
         return sParallelogram(1, s1, s2, b1, b2);
      } else {
         cout << "MpdTpcMiniMC: padANDband's error81" << endl;
         exit(81);
      }
   }

   // d) no crossing (case 4)
   if ((nP1 == 4 && nP2 == 4) || (nM1 == 4 && nM2 == 4)) return 0.;

   cout << "MpdTpcMiniMC: padANDband's error!" << endl;
   exit(90);
}

//______________________________________________________________________
// calculate the rest triangle area of the pad covered by the track band
Double_t TpcMiniMC::sTriangle(Int_t line, Double_t s, Double_t b)
{
   Double_t x0, y0;
   if (TMath::Abs(fEsec.X()) > epsMin && TMath::Abs(fEsec.Y()) > epsMin) {
      if (line == 1) {
         x0 = fR1sec.X() + (b - fR1sec.Y()) * fEsec.X() / fEsec.Y();
         y0 = fR1sec.Y() + (s - fR1sec.X()) * fEsec.Y() / fEsec.X();
      } else {
         x0 = fR2sec.X() + (b - fR2sec.Y()) * fEsec.X() / fEsec.Y();
         y0 = fR2sec.Y() + (s - fR2sec.X()) * fEsec.Y() / fEsec.X();
      }

      // printf("sTriangle: v(%f,%f) line=%d x0=%f y0=%f\n",s,b,line,x0,y0);
      return TMath::Abs((x0 - s) * (y0 - b)) / 2;
   } else {
      printf("sTriangle:: error: there is no crossing with a triangle side\n");
      return 0;
   }
}

//______________________________________________________________________
// calculate the trapezoid area of the pad covered by the track band
Double_t TpcMiniMC::sTrapezoid(Int_t line, Int_t axes, Double_t s, Double_t b1, Double_t b2)
{
   Double_t x1, x2, tng;
   if (axes == 1) { // Y-axes is the trapezoid hight
      tng = fEsec.X() / fEsec.Y();
      if (line == 1) {
         x1 = fR1sec.X() + (b1 - fR1sec.Y()) * tng;
         x2 = fR1sec.X() + (b2 - fR1sec.Y()) * tng;
      } else {
         x1 = fR2sec.X() + (b1 - fR2sec.Y()) * tng;
         x2 = fR2sec.X() + (b2 - fR2sec.Y()) * tng;
      }
   } else {
      tng = fEsec.Y() / fEsec.X();
      if (line == 1) {
         x1 = fR1sec.Y() + (b1 - fR1sec.X()) * tng;
         x2 = fR1sec.Y() + (b2 - fR1sec.X()) * tng;
      } else {
         x1 = fR2sec.Y() + (b1 - fR2sec.X()) * tng;
         x2 = fR2sec.Y() + (b2 - fR2sec.X()) * tng;
      }
   }
   return TMath::Abs((s - 0.5 * (x1 + x2)) * (b2 - b1));
}

//______________________________________________________________________
// calculate the parallelogram area of the pad covered by the track band
Double_t TpcMiniMC::sParallelogram(Int_t axes, Double_t s1, Double_t s2, Double_t b1, Double_t b2)
{
   // if(Nreconstructed==10569) {
   //   printf("sParallelogram: axes=%d s2-s1=%f/%f  b2-b1=%f/%f\n",axes,s2-s1,uEsec.Y(),b2-b1,uEsec.X());
   // }
   if (axes == 1) { // X-axes is the parallelogram hight
      return TMath::Abs(d_track / uEsec.Y() * (s1 - s2));
   } else {
      return TMath::Abs(d_track / uEsec.X() * (b1 - b2));
   }
}

//______________________________________________________________________
// simulate cosmic muon position and direction
void TpcMiniMC::GetCosmicMuon()
{
   Double_t z1, z2, r, theta, phi;
   z1 = -fL_svol;
   z2 = fL_svol;
   if (chambers == 1) {
      z1 = 0;
   } else {
      if (chambers == 2) {
         z2 = 0;
      } else {
         if (chambers != 0) {
            cout << "MpdTpcMiniMC: GetCosmicMuon: no the mode chambers=" << chambers << endl;
            return;
         }
      }
   }
   r   = fr_svol * TMath::Sqrt(fRandom->Rndm());
   phi = TMath::TwoPi() * fRandom->Rndm();
   // random point inside TPC
   r0Beam = TVector3(r * TMath::Cos(phi), r * TMath::Sin(phi), z1 + (z2 - z1) * fRandom->Rndm());
   if (k_par = 2) {
      theta = theta_cosm();
   } else {
      theta = theta_beam();
   }
   phi = TMath::TwoPi() * fRandom->Rndm();
   r   = TMath::Sin(theta),
   // random muon momentum
      pLas = TVector3(r * TMath::Sin(phi), -TMath::Cos(theta), r * TMath::Cos(phi));
   //    pLas=TVector3(0,-1,0);
   //    pLas=TVector3(1,0,0);
   trackPin[0] = r0Beam.X();
   trackPin[1] = r0Beam.Y();
   trackPin[2] = r0Beam.Z();
   trackPin[3] = TMath::ACos(pLas.Z());
   trackPin[4] = TMath::ATan2(pLas.Y(), pLas.X());
   pLas        = (p_min + (p_max - p_min) * fRandom->Rndm()) * pLas;
}

//______________________________________________________________________
// Get the i-th track hit(cluster)
void TpcMiniMC::GetHit(Int_t i, double &x, double &y, double &z, double &w)
{
   x = MpdTpcMiniMC_XHits[i];
   y = MpdTpcMiniMC_YHits[i];
   z = MpdTpcMiniMC_ZHits[i];
   w = MpdTpcMiniMC_WHits[i];
}

//______________________________________________________________________
// Chi2 for the clusters fit by the direct line
double TpcMiniMC_Chi2Line(const double *par)
{
   double chisq = 0;
   double delta;
   Int_t  nhits = MpdTpcMiniMC_Nhits;
   for (int i = 0; i < nhits; i++) {
      double X = MpdTpcMiniMC_XHits[i];
      double Y = MpdTpcMiniMC_YHits[i];
      double Z = MpdTpcMiniMC_ZHits[i];
      double W = MpdTpcMiniMC_WHits[i];

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
   }
   // printf("Chi2Line: chisq=%f  chisq/nh=%f\n",chisq,chisq/MpdTpcMiniMC_Nhits);
   return chisq;
}

//______________________________________________________________________
// find the fit by the direct line
Int_t TpcMiniMC::fitLine()
{
   // Choose method upon creation between:
   // kMigrad, kSimplex, kCombined,
   // kScan, kFumili
   static ROOT::Minuit2::Minuit2Minimizer minFit(ROOT::Minuit2::kMigrad);

   ROOT::Math::Functor fcn(&TpcMiniMC_Chi2Line, 5);
   double              step[5] = {.0001, .0001, .0001, 0.01, 0.01};
   double              limU[5] = {133, 133, 170, pi, pi};
   double              limL[5] = {-133, -133, -170, 0, -pi};
   double              par[5];
   Int_t               imin, imax, iminx, imaxx, iminy, imaxy, iminz, imaxz;
   Double_t            xmn = 200, xmx = -200, ymn = 200, ymx = -200, zmn = 200, zmx = -200, X = 0, Y = 0, Z = 0;
   // printf("fitLine: nHits=%d \n",MpdTpcMiniMC_Nhits);
   for (Int_t i = 0; i < MpdTpcMiniMC_Nhits; i++) {
      if (MpdTpcMiniMC_ZHits[i] > zmx) {
         zmx   = MpdTpcMiniMC_ZHits[i];
         imaxz = i;
      }
      if (MpdTpcMiniMC_ZHits[i] < zmn) {
         zmn   = MpdTpcMiniMC_ZHits[i];
         iminz = i;
      }

      if (MpdTpcMiniMC_XHits[i] > xmx) {
         xmx   = MpdTpcMiniMC_XHits[i];
         imaxx = i;
      }
      if (MpdTpcMiniMC_XHits[i] < xmn) {
         xmn   = MpdTpcMiniMC_XHits[i];
         iminx = i;
      }

      if (MpdTpcMiniMC_YHits[i] > ymx) {
         ymx   = MpdTpcMiniMC_YHits[i];
         imaxy = i;
      }
      if (MpdTpcMiniMC_YHits[i] < ymn) {
         ymn   = MpdTpcMiniMC_YHits[i];
         iminy = i;
      }

      X += MpdTpcMiniMC_XHits[i];
      Y += MpdTpcMiniMC_YHits[i];
      Z += MpdTpcMiniMC_ZHits[i];
   }

   //  for(Int_t i=0;i<MpdTpcMiniMC_Nhits;i++) {
   // printf("fitLine: e=%d i=%d   r(%f,%f,%f) sector=%d\n",Nreconstructed,i,
   // MpdTpcMiniMC_XHits[i],MpdTpcMiniMC_YHits[i],MpdTpcMiniMC_ZHits[i],MpdTpcMiniMC_iSec[i]);
   //}
   Double_t dx = TMath::Abs(MpdTpcMiniMC_XHits[imaxx] - MpdTpcMiniMC_XHits[iminx]);
   Double_t dy = TMath::Abs(MpdTpcMiniMC_YHits[imaxy] - MpdTpcMiniMC_YHits[iminy]);
   Double_t dz = TMath::Abs(MpdTpcMiniMC_ZHits[imaxz] - MpdTpcMiniMC_ZHits[iminz]);
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

   X /= MpdTpcMiniMC_Nhits;
   Y /= MpdTpcMiniMC_Nhits;
   Z /= MpdTpcMiniMC_Nhits;
   TVector3 e =
      TVector3(MpdTpcMiniMC_XHits[imax] - MpdTpcMiniMC_XHits[imin], MpdTpcMiniMC_YHits[imax] - MpdTpcMiniMC_YHits[imin],
               MpdTpcMiniMC_ZHits[imax] - MpdTpcMiniMC_ZHits[imin]);
   if (e.Y() > 0) e = -e;
   e           = e.Unit();
   par[0]      = X;
   par[1]      = Y;
   par[2]      = Z;
   par[3]      = TMath::ACos(e.Z());
   par[4]      = TMath::ATan2(e.Y(), e.X());
   trackPin[5] = 0;
   // printf("par(%f,%f,%f,%f,%f)\n",par[0],par[1],par[2],r2d*par[3],r2d*par[4]);
   for (Int_t i = 0; i < 5; i++) {
      if (i < 3) {
         limL[i] = par[i] - 0.0001;
         limU[i] = par[i] + 0.0001;
      } else {
         limL[i] = par[i] - 0.1 * TMath::Abs(par[i]);
         limU[i] = par[i] + 0.1 * TMath::Abs(par[i]);
      }
   };

   minFit.SetFunction(fcn);

   // Set the free variables to be minimized!
   minFit.SetLimitedVariable(0, "x0", par[0], step[0], limL[0], limU[0]);
   minFit.SetLimitedVariable(1, "y0", par[1], step[1], limL[1], limU[1]);
   minFit.SetLimitedVariable(2, "z0", par[2], step[2], limL[2], limU[2]);
   minFit.SetLimitedVariable(3, "theta", par[4], step[3], limL[3], limU[3]);
   minFit.SetLimitedVariable(4, "phi", par[4], step[4], limL[4], limU[4]);

   Bool_t res = minFit.Minimize();
   // printf("fitLine:  event=%6d yes=%d  E/n=%9.5f\n",Nreconstructed,res,minFit.MinValue()/MpdTpcMiniMC_Nhits);
   // if(!res) minFit.PrintResults();
   //     if(!res) {minFit.PrintResults(); }
   // printf("fitLine: minfit/nhits=%f EvalpHit=%f \n",minFit.MinValue()/MpdTpcMiniMC_Nhits,EvalpHit);
   if (minFit.MinValue() / MpdTpcMiniMC_Nhits > EvalpHit) return 1;
   //      minFit.PrintResults();
   const double *xs = minFit.X();
   trackPout[0]     = xs[0];
   trackPout[1]     = xs[1];
   trackPout[2]     = xs[2];
   if (TMath::Abs(trackPin[3] - xs[3]) > 1)
      trackPout[3] = pi - xs[3];
   else
      trackPout[3] = xs[3];
   if ((xs[4] - trackPin[4]) > 1)
      trackPout[4] = xs[4] - pi;
   else {
      if ((trackPin[4] - xs[4]) < -1)
         trackPout[4] = xs[4] + pi;
      else
         trackPout[4] = xs[4];
   }
   //      trackPout[3]=xs[3];
   //      trackPout[4]=xs[4];
   trackPout[5] = minFit.NCalls();
   // printf("Par: 0(%9.4f %9.4f) 1(%9.4f %9.4f) 2(%9.4f %9.4f) 3(%9.4f %9.4f) 4(%9.4f %9.4f)\n",
   // trackPin[0],trackPout[0],trackPin[1],trackPout[1],trackPin[2],trackPout[2],r2d*trackPin[3],r2d*trackPout[3],r2d*trackPin[4],r2d*trackPout[4]);
   for (Int_t k = 0; k < 6; k++) trackDiff[k] = trackPout[k] - trackPin[k];
   return 0;
   //    }
   //    else return 1;
}

//______________________________________________________________________
// put chi2 values
void TpcMiniMC::putChi2L()
{

   Double_t chisq2 = 0;
   TVector3 e2(TMath::Sin(trackPin[3]) * TMath::Cos(trackPin[4]), TMath::Sin(trackPin[3]) * TMath::Sin(trackPin[4]),
               TMath::Cos(trackPin[3]));
   TVector3 x02(trackPin[0], trackPin[1], trackPin[2]);
   //    x02=TVector3(trackPout[0],trackPout[1],trackPout[2]);
   TVector3 dx2;
   e2 = 1 / pLas.Mag() * pLas;

   Double_t chisq = 0;
   Double_t delta;
   TVector3 x1, dx;
   Int_t    nhits = MpdTpcMiniMC_Nhits;
   Double_t ex    = TMath::Sin(trackPout[3]) * TMath::Cos(trackPout[4]);
   Double_t ey    = TMath::Sin(trackPout[3]) * TMath::Sin(trackPout[4]);
   Double_t ez    = TMath::Cos(trackPout[3]);
   TVector3 e     = TVector3(ex, ey, ez);
   TVector3 x0    = TVector3(trackPout[0], trackPout[1], trackPout[2]);
   // printf("iChi2: x0(%f,%f,%f) e(%f,%f,%f)\n",x0.X(),x0.Y(),x0.Z(),e.X(),e.Y(),e.Z());
   // printf("oChi2: x02(%f,%f,%f) e2(%f,%f,%f)\n",x02.X(),x02.Y(),x02.Z(),e2.X(),e2.Y(),e2.Z());
   for (int i = 0; i < nhits; i++) {
      Double_t X = MpdTpcMiniMC_XHits[i];
      Double_t Y = MpdTpcMiniMC_YHits[i];
      Double_t Z = MpdTpcMiniMC_ZHits[i];
      Double_t W = MpdTpcMiniMC_WHits[i];
      x1         = TVector3(MpdTpcMiniMC_XHits[i], MpdTpcMiniMC_YHits[i], MpdTpcMiniMC_ZHits[i]);
      // printf("putChi2:  (x0-x1)=(%f,%f,%f)\n",x0.X()-x1.X(),x0.Y()-x1.Y(),x0.Z()-x1.Z());
      dx                    = e.Cross(x0 - x1);
      delta                 = dx.Mag() / W;
      MpdTpcMiniMC_chi2s[i] = delta * delta;
      chisq += delta * delta;
      // printf("putChi2: MpdTpcMiniMC_iSec[%d]=%d delta^2=%f chi2=%f\n",i,MpdTpcMiniMC_iSec[i],delta*delta,chisq);
      // printf("nputChi2: (x02-x1)=(%f,%f,%f)\n",x02.X()-x1.X(),x02.Y()-x1.Y(),x02.Z()-x1.Z());
      dx2   = e2.Cross(x02 - x1);
      delta = dx2.Mag() / W;
      chisq2 += delta * delta;
      // printf("putChi2: MpdTpcMiniMC_iSec[%d] delta^2=%f chi2=%f\n\n",i,delta*delta,chisq2);

      // printf("putChi2: MpdTpcMiniMC_iSec[%d]=%d\n",i,MpdTpcMiniMC_iSec[i]);
      fNhitSec[MpdTpcMiniMC_iSec[i]]++;
      fSecChi2[MpdTpcMiniMC_iSec[i]] += delta * delta;
   }
   Chi2Fit = chisq / MpdTpcMiniMC_Nhits;
   // printf("putChi2: Chi2Fit=(%f,%f) tPout3,4=(%f,%f) tPin3,4=(%f,%f)\n",
   // Chi2Fit,chisq2/MpdTpcMiniMC_Nhits,r2d*trackPout[3],r2d*trackPout[4],r2d*trackPin[3],r2d*trackPin[4]);
}

//  *************************************************************************
//  ************************** magnetic field part **************************
//  *************************************************************************
//______________________________________________________________________
// Simulate muon parameters at a magnetic field
void TpcMiniMC::GetPointMuon()
{
   Double_t theta1, theta2, r, theta, phi, pT, rc2B, rcB, phi_s, t, t1, t2, arg;
   q_muon = fRandom->Rndm() > 0.5 ? 1 : -1;
   q_muon = -q_muon;
   q_sign = q_muon / TMath::Abs(q_muon);
   // rotation sign for helix
   signRotB = (q_sign > 0) ? -1 : 1;
   // muon momentum
   theta1 = minMuonTheta;
   theta2 = TMath::Pi() - minMuonTheta;
   if (chambers == 1) {
      theta2 = TMath::Pi() / 2;
   } else {
      if (chambers == 2) {
         theta1 = TMath::Pi() / 2;
      } else {
         if (chambers != 0) {
            cout << "MpdTpcMiniMC: GetCosmicMuon: no the mode chambers=" << chambers << endl;
            return;
         }
      }
   }
   r0BeamB = G2B * r0Beam; // production point in BCS
   theta   = theta1 + (theta2 - theta1) * fRandom->Rndm();
   //    theta=40.914114*d2r;
   phi = TMath::TwoPi() * fRandom->Rndm();
   //    phi=75.085936*d2r;
   // printf("  %f  %f\n",   theta*r2d , phi*r2d);
   pT   = p_min + (p_max - p_min) * fRandom->Rndm();
   pLas = TVector3(TMath::Sin(theta) * TMath::Cos(phi), TMath::Sin(theta) * TMath::Sin(phi), TMath::Cos(theta));
   pLas = (TMath::Abs(theta) > epsMin) ? pT / TMath::Sin(theta) * pLas : 1e10 * pLas;
   Double_t pMag = pLas.Mag();
   pMuonB        = G2B * pLas;
   pTMuonB       = TMath::Sqrt(pMuonB.X() * pMuonB.X() + pMuonB.Y() * pMuonB.Y());
   Double_t bMag = bMf.Mag();
   // Helix parameters
   rH       = 333.56 * pTMuonB / (TMath::Abs(q_muon) * bMag); // base helix radius
   rHa[0]   = rH - dh_track;                                  // inner helix radius
   rHa[1]   = rH + dh_track;                                  // outer helix radius
   rHa[2]   = rH;                                             // outer helix radius
   xcB      = r0BeamB.X() + rH * pMuonB.Y() / (q_muon * pTMuonB);
   ycB      = r0BeamB.Y() - rH * pMuonB.X() / (q_muon * pTMuonB);
   rc2B     = xcB * xcB + ycB * ycB;
   rcB      = TMath::Sqrt(rc2B);
   rc0_B    = TVector3(xcB, ycB, 0); // point on the cilinder axis in BSC
   rc0_G    = B2G * rc0_B;           // point on the cilinder axis in GSC
   phi_0B   = TMath::ATan2(-ycB, -xcB);
   helix_Pz = rH * pMuonB.Z() / pTMuonB;
   // printf("(parHelix(%6.1f,%6.1f,%f,%6.1f,%f)\n",xcB,ycB,rH,phi_0B*r2d,helix_Pz);
   //  muon end point in the chamber
   //  cross point with chamber barrel
   t   = 100;
   arg = fr_svol * fr_svol / (2 * rc2B) - 1;
   if (TMath::Abs(arg) < 1) { // helix crosses the chamber cilinder
      phi_s       = TMath::ATan2(xcB, ycB);
      Double_t a1 = TMath::ASin(arg);
      Double_t a2 = (arg > 0) ? pi - a1 : -pi - a1;
      t1          = signRotB * (a1 - phi_0B - phi_s);
      t2          = signRotB * (a2 - phi_0B - phi_s);
      if (t1 > twopi)
         t1 -= twopi;
      else {
         if (t1 < -twopi) t1 += twopi;
      };
      if (t1 < 0) t1 += twopi;
      if (t2 > twopi)
         t2 -= twopi;
      else {
         if (t2 < -twopi) t2 += twopi;
      };
      if (t2 < 0) t2 += twopi;
      if (t1 < t2)
         t = t1;
      else
         t = t2;
   }
   if (pMuonB.Z() == 0) {
      if (t == 100)
         t_end = twopi - phi_sec;
      else
         t_end = t;
   } else {
      t1 = pMuonB.Z() / TMath::Abs(pMuonB.Z());
      t2 = (t1 * fL_svol - r0Beam.Z()) / (TMath::Abs(bMf.Z() / bMf.Mag()) * helix_Pz);
      if (t2 < twopi || t < twopi) {
         if (t2 < t)
            t_end = t2;
         else
            t_end = t;
      } else {
         //        t_end=twopi-phi_sec;
         t_end = twopi;
      }
   }
   // cout<<"("<<t<<","<<t2<<")  t_end="<<t_end<<endl;
   trackPin[0] = xcB;
   trackPin[1] = ycB;
   trackPin[2] = r0BeamB.Z();
   trackPin[3] = q_muon * rH;
   trackPin[4] = phi_0B;
   trackPin[5] = helix_Pz;
}

//______________________________________________________________________
// simulate pads signals for a muon in the magnetic field
Int_t TpcMiniMC::passMuoninMF()
{
   // printf("MpdTpcMiniMC::passMuoninMF()  t_end=%f================\n",t_end);
   Bool_t   InPad;
   Int_t    nclust, irow, irow0, ipad, ipad0, ip1, ip2, nPads, iSector0, sectorI, sec1;
   Double_t t, t0, ts, t0_end, dt, step = 0.1, tc, xt, xt0, xc, x1, x2, yt, yt0, yc, y1, y2, hPad, wPad;
   MpdTpcMiniMC_Nhits = 0;
   dt                 = step / rH; // step along phi
   t = ts   = 0;
   iSector0 = -20;
   t0_end   = t_end;
   InPad    = false;   // particle did not cross any sensitive sector area
   while (t < t_end) { // loop along full helix
                       // printf("while(t<=t_end)   t=%f  t_end=%f tc=%f\n",t*r2d,t_end*r2d,tc*r2d);
      for (Int_t i = 0; i < mxp; i++)
         for (Int_t j = 0; j < mxr; j++)
            for (Int_t k = 0; k < 3; k++) indexPad[i][j][k] = 0;
      sec1 = 0;
      //    wasInPad=false;
      for (Int_t l = 2; l >= 0; l--) { // 3 trips along the helix
                                       // l=2,1,0 own particle, inner and outer circles
         if (l == 2) {
            t0     = ts;
            t0_end = t_end;
         } else {
            if (sec1 == 0)
               t0_end = t_end;
            else
               t0_end = ts;
         }
         t = t0; // start point
         // printf("for  l=%d ===== t=%f t0_end=%f =========\n",l,t*r2d,t0_end*r2d);
         iSector  = 25;
         iSector0 = coordinates(t, rHa[l]);
         ipad     = -1;
         irow     = -1;
         while (t <= t0_end) {
            rLt0 = rLt;
            t += dt;
            iSector = coordinates(t, rHa[l]);
            // printf("while(t=%f== l=%d ==  sec=%d sec0=%d \n",t*r2d,l,iSector,iSector0);
            if (iSector == iSector0 || l != 2) {
               xt0   = rLt0.X();
               yt0   = rLt0.Y();
               xt    = rLt.X();
               yt    = rLt.Y();
               ipad0 = ipad;
               irow0 = irow;
               if (fTpcSecGeo->IsItPad(xt, yt, irow, ipad, x1, x2, y1, y2)) {
                  x1 = xminPad[ipad][irow];
                  y1 = yminPad[ipad][irow];
                  x2 = xmaxPad[ipad][irow];
                  y2 = ymaxPad[ipad][irow];
                  // the point is inside the pad(irow,ipad)
                  InPad = true;
                  if (l == 2) sectorI = iSector;
                  if (ipad != ipad0 || irow != irow0) { // just a NEXT pad
                                                        // if(l==2)
                     // printf(" IsItPad l=%d t=%f sec=%d i,p(%d,%d) rBt(%f,%f,%f) rLt(%f,%f)\n",
                     // l,t*r2d,iSector,ipad,irow,rBt.X(),rBt.Y(),rBt.Z(),xt,yt);
                     //  check the bottom border crossing
                     tc = (y1 - yt0) / (yt - yt0); // [rt0,rt] line parameter t at y=y1
                     xc = xt0 + (xt - xt0) * tc;
                     if ((rLt0.Y() <= y1 && y1 <= rLt.Y() || rLt.Y() <= y1 && y1 <= rLt0.Y()) && x1 <= xc &&
                         xc <= x2) { // track crosses the lower bound
                        if (l < 2) {
                           xc1[ipad][irow][l] = xc;
                           yc1[ipad][irow][l] = y1;
                        }
                        indexPad[ipad][irow][l] = 1;          // bottom border cross index
                        if (irow0 >= 0) {                     // the previous point was inside other pad
                           if (xc <= xminPad[ipad0][irow0]) { // came from r-b pad
                              if (l < 2) {
                                 tc                   = (x2 - xt0) / (xt - xt0);
                                 yc                   = yt0 + (yt - yt0) * tc;
                                 xc2[ipad0][irow0][l] = x2;
                                 yc2[ipad0][irow0][l] = yc;
                              }
                              indexPad[ipad0][irow0][l] += 1000; // left border cross index
                           } else {
                              if (xc >= xmaxPad[ipad0][irow0]) { // came from l-b pad
                                 if (l < 2) {
                                    tc                   = (x1 - xt0) / (xt - xt0);
                                    yc                   = yt0 + (yt - yt0) * tc;
                                    xc2[ipad0][irow0][l] = x1;
                                    yc2[ipad0][irow0][l] = yc;
                                 }
                                 indexPad[ipad0][irow0][l] += 10; // right border cross index
                              } else {                            // came from bottom pad
                                 if (l < 2) {
                                    xc2[ipad0][irow0][l] = xc;
                                    yc2[ipad0][irow0][l] = y1;
                                 }
                                 indexPad[ipad0][irow0][l] += 100; // top border cross index
                              }
                           }
                        }
                     } else {
                        // check the top border crossing
                        tc = (y2 - yt0) / (yt - yt0); // [rt0,rt] lleftine parameter t at y=y2
                        xc = xt0 + (xt - xt0) * tc;
                        if ((rLt0.Y() <= y2 && y2 <= rLt.Y() || rLt.Y() <= y2 && y2 <= rLt0.Y()) && x1 <= xc &&
                            xc <= x2) { // track crosses the top bound
                           if (l < 2) {
                              xc1[ipad][irow][l] = xc;
                              yc1[ipad][irow][l] = y2;
                           }
                           indexPad[ipad][irow][l] = 100;        // top border cross index
                           if (irow0 >= 0) {                     // the previous point was inside other pad
                              if (xc <= xminPad[ipad0][irow0]) { // came from r-b pad
                                 if (l < 2) {
                                    tc                   = (x2 - xt0) / (xt - xt0);
                                    yc                   = yt0 + (yt - yt0) * tc;
                                    xc2[ipad0][irow0][l] = x2;
                                    yc2[ipad0][irow0][l] = yc;
                                 }
                                 indexPad[ipad0][irow0][l] += 1000; // l-border cross index
                              } else {
                                 if (xc >= xmaxPad[ipad0][irow0]) { // came from l-b pad
                                    if (l < 2) {
                                       tc                   = (x1 - xt0) / (xt - xt0);
                                       yc                   = yt0 + (yt - yt0) * tc;
                                       xc2[ipad0][irow0][l] = x1;
                                       yc2[ipad0][irow0][l] = yc;
                                    }
                                    indexPad[ipad0][irow0][l] += 10; // r-border cross index
                                 } else {                            // came from bottom pad
                                    if (l < 2) {
                                       xc2[ipad0][irow0][l] = xc;
                                       yc2[ipad0][irow0][l] = y2;
                                    }
                                    indexPad[ipad0][irow0][l] += 1; // b-border cross index
                                 }
                              }
                           }
                        } else {
                           // check the left border crossing
                           tc = (x1 - xt0) / (xt - xt0); // [rt0,rt] line parameter t at x=x1
                           yc = yt0 + (yt - yt0) * tc;
                           if ((rLt0.X() <= x1 && x1 <= rLt.X() || rLt.X() <= x1 && x1 <= rLt0.X()) && y1 <= yc &&
                               yc <= y2) { // track crosses the left border
                              if (l < 2) {
                                 xc1[ipad][irow][l] = x1;
                                 yc1[ipad][irow][l] = yc;
                              }
                              indexPad[ipad][irow][l] = 1000;       // left border cross index
                              if (irow0 >= 0) {                     // the previous point was inside other pad
                                 if (yc <= yminPad[ipad0][irow0]) { // came from l-t pad
                                    if (l < 2) {
                                       tc                   = (y2 - yt0) / (yt - yt0);
                                       xc                   = xt0 + (xt - xt0) * tc;
                                       xc2[ipad0][irow0][l] = xc;
                                       yc2[ipad0][irow0][l] = y2;
                                    }
                                    indexPad[ipad0][irow0][l] += 1; // b-border cross index
                                 } else {
                                    if (yc >= ymaxPad[ipad0][irow0]) { // came from l-b pad
                                       if (l < 2) {
                                          tc                   = (y1 - yt0) / (yt - yt0);
                                          xc                   = xt0 + (xt - xt0) * tc;
                                          xc2[ipad0][irow0][l] = xc;
                                          yc2[ipad0][irow0][l] = y1;
                                       }
                                       indexPad[ipad0][irow0][l] += 100; // t-border cross index
                                    } else {                             // came from upper pad
                                       if (l < 2) {
                                          xc2[ipad0][irow0][l] = x1;
                                          yc2[ipad0][irow0][l] = yc;
                                       }
                                       indexPad[ipad0][irow0][l] += 10; // r-border cross index
                                    }
                                 }
                              }
                           } else {
                              // check the right border crossing
                              tc = (x2 - xt0) / (xt - xt0); // [rt0,rt] line parameter t at x=x2
                              yc = yt0 + (yt - yt0) * tc;
                              if ((rLt0.X() <= x2 && x2 <= rLt.X() || rLt.X() <= x2 && x2 <= rLt0.X()) && y1 <= yc &&
                                  yc <= y2) { // track crosses the right border
                                 if (l < 2) {
                                    xc1[ipad][irow][l] = x2;
                                    yc1[ipad][irow][l] = yc;
                                 }
                                 indexPad[ipad][irow][l] = 10;         // right border cross index
                                 if (irow0 >= 0) {                     // the previous point was inside other pad
                                    if (yc <= yminPad[ipad0][irow0]) { // came from r-b pad
                                       if (l < 2) {
                                          tc                   = (y2 - yt0) / (yt - yt0);
                                          xc                   = xt0 + (xt - xt0) * tc;
                                          xc2[ipad0][irow0][l] = xc;
                                          yc2[ipad0][irow0][l] = y2;
                                       }
                                       indexPad[ipad0][irow0][l] += 1; // b-border cross index
                                    } else {
                                       if (yc >= ymaxPad[ipad0][irow0]) { // came from l-b pad
                                          if (l < 2) {
                                             tc                   = (y1 - yt0) / (yt - yt0);
                                             xc                   = xt0 + (xt - xt0) * tc;
                                             xc2[ipad0][irow0][l] = xc;
                                             yc2[ipad0][irow0][l] = y1;
                                          }
                                          indexPad[ipad0][irow0][l] += 100; // top border index
                                       } else {                             // came from upper pad
                                          if (l < 2) {
                                             xc2[ipad0][irow0][l] = x1;
                                             yc2[ipad0][irow0][l] = yc;
                                          }
                                          indexPad[ipad0][irow0][l] += 1000; // left border index
                                       }
                                    }
                                 }
                              } else {
                                 printf("error at the check of the border crossing\n");
                              }
                              //                    ipad0=ipad;
                              //                    irow0=irow;
                           } // else of r-b crossing
                        }    // else of l-b crossing
                     }       // else of u-b crossing
                  }          // check of ipad0!=ipad
               }             // check of IsItPad
               else {        // track is out of pad
                  // printf("0000000 l=%d 00000   t=%f   rt0(%f,%f) t0_end=%f sec=%d\n",
                  //     l,t*r2d,xt0,yt0,t0_end*r2d,iSector);
                  //             if(InPad) InPad=false;
                  if (ipad0 > 0) { // exit from the pad(ipad0,irow0)
                     // find the border of the exit
                     // x1,y1,x2,y2 for pad(ipad0,irow0)
                     x1 = xminPad[ipad0][irow0];
                     y1 = yminPad[ipad0][irow0];
                     x2 = xmaxPad[ipad0][irow0];
                     y2 = ymaxPad[ipad0][irow0];
                     // printf("p,r(%d,%d) rt0(%f,%f) rt(%f,%f) r1(%f,%f) r2(%f,%f)\n",
                     // ipad0,irow0,xt0,yt0,xt,yt,x1,y1,x2,y2);
                     //  check the bottom border crossing
                     tc = (y1 - yt0) / (yt - yt0); // [rt0,rt] line parameter t at y=y1
                     xc = xt0 + (xt - xt0) * tc;
                     if ((rLt0.Y() <= y1 && y1 <= rLt.Y() || rLt.Y() <= y1 && y1 <= rLt0.Y()) && x1 <= xc &&
                         xc <= x2) { // track crosses the lower bound
                        if (l < 2) {
                           xc2[ipad0][irow0][l] = xc;
                           yc2[ipad0][irow0][l] = y1;
                        }
                        indexPad[ipad0][irow0][l] += 1; // bottom border cross index
                     } else {
                        // check the top border crossing
                        tc = (y2 - yt0) / (yt - yt0); // [rt0,rt]     line parameter t at y=y2
                        xc = xt0 + (xt - xt0) * tc;
                        if ((rLt0.Y() <= y2 && y2 <= rLt.Y() || rLt.Y() <= y2 && y2 <= rLt0.Y()) && x1 <= xc &&
                            xc <= x2) { // track crosses the top bound
                           if (l < 2) {
                              xc2[ipad0][irow0][l] = xc;
                              yc2[ipad0][irow0][l] = y2;
                           }
                           indexPad[ipad0][irow0][l] += 100; // top border cross index
                        } else {
                           // check the left border crossing
                           tc = (x1 - xt0) / (xt - xt0); // [rt0,rt] line parameter t at x=x1
                           yc = yt0 + (yt - yt0) * tc;
                           if ((rLt0.X() <= x1 && x1 <= rLt.X() || rLt.X() <= x1 && x1 <= rLt0.X()) && y1 <= yc &&
                               yc <= y2) { // track crosses the left border
                              if (l < 2) {
                                 xc2[ipad0][irow0][l] = x1;
                                 yc2[ipad0][irow0][l] = yc;
                              }
                              indexPad[ipad0][irow0][l] += 1000; // l-border cross index
                           } else {
                              // check the right border crossing
                              tc = (x2 - xt0) / (xt - xt0); // [rt0,rt] line parameter t at x=x2
                              yc = yt0 + (yt - yt0) * tc;
                              if ((rLt0.X() <= x2 && x2 <= rLt.X() || rLt.X() <= x2 && x2 <= rLt0.X()) && y1 <= yc &&
                                  yc <= y2) { // track crosses the right border
                                 if (l < 2) {
                                    xc1[ipad0][irow0][l] = x2;
                                    yc1[ipad0][irow0][l] = yc;
                                 }
                                 indexPad[ipad0][irow0][l] = +10; // r-border cross index
                              } else
                                 printf("error 2 at exit from pad in MuonPass\n");
                           }
                        }
                     }
                     ipad = -1;
                     irow = -1;
                  }
                  continue;
               }
            }      // check of iSector==iSector0
            else { // new sector
                   //          sectorI=iSector0;
                   //          iSector0=iSector;
               if (l == 2) {
                  if (InPad) {
                     sectorI = iSector0;
                     ts      = t;
                     sec1    = 1;
                     break;
                  } else {
                     t0       = t;
                     iSector0 = iSector;
                     continue;
                  }
               } else {
                  // printf("l=%d t=%f sec1=%d t0=%f t0_end=%f iSector=%d iSector0=%d\n",
                  // l,t*r2d,sec1,t0*r2d,t0_end*r2d,iSector,iSector0);
                  //             if(iSector!=sec1) continue; else break;
                  iSector0 = iSector;
                  continue;
               }
            }
         } // end of loop while(t<=t0_end)
           //      if(l==2&&t>t_end) sectorI=iSector;
           // cout<<"======= l="<<l<<" ==========2=================t="<<t*r2d<<endl;
      }    // end of loop for(Int_t l=2;l>=0;l--)
           // printf(" t=%f t_end=%f ts=%f wasInPad=%d iSector=%d sectorI=%d =====================\n",
      // t*r2d,t_end*r2d,ts*r2d,wasInPad,iSector,sectorI);
      //     if(t<t_end) iSector=sectorI;
      //     if(ts>t_end) iSector=sectorI;
      //     iSector=sectorI;
      //       if(sec1!=0)
      iSector = sectorI;
      fTpcSecGeo->SectorTransformation(iSector, R0shift, loc2glo);
      e_3 = TVector3(loc2glo.XZ(), loc2glo.YZ(), loc2glo.ZZ());
      fTpcSecGeo->SectorBackTransformation(iSector, G0shift, glo2loc);

      /*for(Int_t k=0;k<3;k++) {
        printf("========================================= %d =================\n",k);
        for(Int_t j=0;j<mxr;j++) {
          for(Int_t i=0;i<mxp;i++)
            if(indexPad[i][j][k]!=0) printf("X"); else printf("0");
          cout<<endl;
        }
      }*/
      //  printf("========================================= %d =================\n",iSector);
      CheckPadsinMF();
      // printf("end of loop for t=%f sec=%d hits=%d\n",t*r2d,iSector,MpdTpcMiniMC_Nhits);
   } // end of loop while(t<=t0)
   if (MpdTpcMiniMC_Nhits >= fminHits)
      return 1;
   else
      return 0;
}

//______________________________________________________________________
// find all clusters in a sector from a particle in the magnetic field
void TpcMiniMC::CheckPadsinMF()
{
   Int_t    nclust = 0, irow, ipad, ip1, ip2, nPads, nRows;
   Double_t hPad, wPad, sw1, sw2, x1, y1, x2, y2;
   //   ***************************************************
   //   ********       check pads in each row      ********
   //   ***************************************************

   wPad  = fWpad;
   nRows = fTpcSecGeo->NofRows();
   // printf("==================== Sector=%d ======================\n",iSector);
   for (irow = 0; irow < nRows; irow++) { //  find clusters in the row
                                          // printf("passMuoninMF: *************** irow=%d *****************\n",irow);
      hPad  = fTpcSecGeo->PadHeight(irow);
      nPads = fTpcSecGeo->NofPadsInRow(irow);
      ip1   = 0;
      ip2   = ip1 + nPads;
      for (ipad = ip1; ipad < ip2; ipad++) {
         x1 = xminPad[ipad][irow];
         y1 = yminPad[ipad][irow];
         x2 = xmaxPad[ipad][irow];
         y2 = ymaxPad[ipad][irow];

         // if(indexPad[ipad][irow][0]!=0||indexPad[ipad][irow][1]!=0
         //||indexPad[ipad][irow][2]!=0)
         // printf("CheckPadsinMF: sec=%2d indexPad[%2d,%2d](0,1,2)=(%4d %4d %4d)\n",iSector,
         // ipad,irow,indexPad[ipad][irow][0],indexPad[ipad][irow][1],indexPad[ipad][irow][2]);

         if (indexPad[ipad][irow][0] == 0 && indexPad[ipad][irow][1] == 0) {
            if (indexPad[ipad][irow][2] != 0) { // case a) of padANDband
               SPads[ipad] = wPad * hPad;
            } else {
               SPads[ipad] = 0;
            }
            continue;
         }

         //                     base 1100
         // case b) or d)
         if (indexPad[ipad][irow][0] == 1100 && indexPad[ipad][irow][1] == 0) {
            if (xc1[ipad][irow][0] < xc2[ipad][irow][0])
               sw1 = 0.5 * (y2 - yc1[ipad][irow][0]) * (xc2[ipad][irow][0] - x1);
            else
               sw1 = 0.5 * (y2 - yc2[ipad][irow][0]) * (xc1[ipad][irow][0] - x1);
            if (ipad < ip2 - 1)
               if (indexPad[ipad + 1][irow][1] != 0)
                  SPads[ipad] = wPad * hPad - sw1; // b)
               else
                  SPads[ipad] = sw1; // d)
            else if (indexPad[ipad - 1][irow][1] != 0)
               SPads[ipad] = sw1; // d)
            else
               SPads[ipad] = wPad * hPad - sw1; // b)
            continue;
         }
         if (indexPad[ipad][irow][1] == 1100 && indexPad[ipad][irow][0] == 0) {
            if (xc1[ipad][irow][1] < xc2[ipad][irow][1])
               sw1 = 0.5 * (y2 - yc1[ipad][irow][1]) * (xc2[ipad][irow][1] - x1);
            else
               sw1 = 0.5 * (y2 - yc2[ipad][irow][1]) * (xc1[ipad][irow][1] - x1);
            if (ipad < ip2 - 1)
               if (indexPad[ipad + 1][irow][0] != 0)
                  SPads[ipad] = wPad * hPad - sw1; // b)
               else
                  SPads[ipad] = sw1; // d)
            else if (indexPad[ipad - 1][irow][0] != 0)
               SPads[ipad] = sw1; // d)
            else
               SPads[ipad] = wPad * hPad - sw1; // b)
            continue;
         }
         // case c)
         if (indexPad[ipad][irow][0] == 101 && indexPad[ipad][irow][1] == 0) {
            if (ipad < ip2 - 1)
               if (indexPad[ipad + 1][irow][1] != 0)
                  SPads[ipad] = (x2 - 0.5 * (xc1[ipad][irow][0] + xc2[ipad][irow][0])) * hPad;
               else
                  SPads[ipad] = (0.5 * (xc1[ipad][irow][0] + xc2[ipad][irow][0]) - x1) * hPad;
            else if (indexPad[ipad - 1][irow][1] != 0)
               SPads[ipad] = (0.5 * (xc1[ipad][irow][0] + xc2[ipad][irow][0]) - x1) * hPad;
            else
               SPads[ipad] = (x2 - 0.5 * (xc1[ipad][irow][0] + xc2[ipad][irow][0])) * hPad;
            continue;
         }
         if (indexPad[ipad][irow][1] == 101 && indexPad[ipad][irow][0] == 0) {
            if (ipad < ip2 - 1)
               if (indexPad[ipad + 1][irow][0] != 0)
                  SPads[ipad] = (x2 - 0.5 * (xc1[ipad][irow][1] + xc2[ipad][irow][1])) * hPad;
               else
                  SPads[ipad] = (0.5 * (xc1[ipad][irow][1] + xc2[ipad][irow][1]) - x1) * hPad;
            else if (indexPad[ipad - 1][irow][0] != 0)
               SPads[ipad] = (0.5 * (xc1[ipad][irow][1] + xc2[ipad][irow][1]) - x1) * hPad;
            else
               SPads[ipad] = (x2 - 0.5 * (xc1[ipad][irow][1] + xc2[ipad][irow][1])) * hPad;
            continue;
         }
         // case e)
         if (indexPad[ipad][irow][0] == 1100 && indexPad[ipad][irow][1] == 11) {
            if (xc1[ipad][irow][0] == x1) {
               sw1 = 0.5 * (y2 - yc1[ipad][irow][0]) * (xc2[ipad][irow][0] - x1);
               if (yc1[ipad][irow][1] == y1)
                  sw2 = 0.5 * (x2 - xc1[ipad][irow][1]) * (yc2[ipad][irow][1] - y1);
               else
                  sw2 = 0.5 * (x2 - xc2[ipad][irow][1]) * (yc1[ipad][irow][1] - y1);
            } else {
               sw1 = 0.5 * (y2 - yc2[ipad][irow][0]) * (xc1[ipad][irow][0] - x1);
               if (yc1[ipad][irow][1] == y1)
                  sw2 = 0.5 * (x2 - xc1[ipad][irow][1]) * (yc2[ipad][irow][1] - y1);
               else
                  sw2 = 0.5 * (x2 - xc2[ipad][irow][1]) * (yc1[ipad][irow][1] - y1);
            }
            SPads[ipad] = wPad * hPad - sw1 - sw2;
            continue;
         }
         if (indexPad[ipad][irow][1] == 1100 && indexPad[ipad][irow][0] == 11) {
            if (xc1[ipad][irow][1] == x1) {
               sw1 = 0.5 * (y2 - yc1[ipad][irow][1]) * (xc2[ipad][irow][1] - x1);
               if (yc1[ipad][irow][0] == y1)
                  sw2 = 0.5 * (x2 - xc1[ipad][irow][0]) * (yc2[ipad][irow][0] - y1);
               else
                  sw2 = 0.5 * (x2 - xc2[ipad][irow][0]) * (yc1[ipad][irow][0] - y1);
            } else {
               sw1 = 0.5 * (y2 - yc2[ipad][irow][1]) * (xc1[ipad][irow][1] - x1);
               if (yc1[ipad][irow][0] == y1)
                  sw2 = 0.5 * (x2 - xc1[ipad][irow][0]) * (yc2[ipad][irow][0] - y1);
               else
                  sw2 = 0.5 * (x2 - xc2[ipad][irow][0]) * (yc1[ipad][irow][0] - y1);
            }
            SPads[ipad] = wPad * hPad - sw1 - sw2;
            continue;
         }
         // case f)
         if (indexPad[ipad][irow][0] == 1100 && indexPad[ipad][irow][1] == 101) {
            if (xc1[ipad][irow][0] == x1)
               sw1 = 0.5 * (y2 - yc1[ipad][irow][0]) * (xc2[ipad][irow][0] - x1);
            else
               sw1 = 0.5 * (y2 - yc2[ipad][irow][0]) * (xc1[ipad][irow][0] - x1);
            sw2         = (x2 - 0.5 * (xc1[ipad][irow][1] + xc2[ipad][irow][1])) * hPad;
            SPads[ipad] = wPad * hPad - sw1 - sw2;
            continue;
         }
         if (indexPad[ipad][irow][1] == 1100 && indexPad[ipad][irow][0] == 101) {
            if (xc1[ipad][irow][1] == x1)
               sw1 = 0.5 * (y2 - yc1[ipad][irow][1]) * (xc2[ipad][irow][1] - x1);
            else
               sw1 = 0.5 * (y2 - yc2[ipad][irow][1]) * (xc1[ipad][irow][1] - x1);
            sw2         = (x2 - 0.5 * (xc1[ipad][irow][0] + xc2[ipad][irow][0])) * hPad;
            SPads[ipad] = wPad * hPad - sw1 - sw2;
            continue;
         }
         if (indexPad[ipad][irow][0] == 1100 && indexPad[ipad][irow][1] == 1010) {
            if (xc1[ipad][irow][0] == x1)
               sw1 = 0.5 * (y2 - yc1[ipad][irow][0]) * (xc2[ipad][irow][0] - x1);
            else
               sw1 = 0.5 * (y2 - yc2[ipad][irow][0]) * (xc1[ipad][irow][0] - x1);
            sw2         = (0.5 * (yc1[ipad][irow][1] + yc2[ipad][irow][1]) - y1) * wPad;
            SPads[ipad] = wPad * hPad - sw1 - sw2;
            continue;
         }
         if (indexPad[ipad][irow][1] == 1100 && indexPad[ipad][irow][0] == 1010) {
            if (xc1[ipad][irow][1] == x1)
               sw1 = 0.5 * (y2 - yc1[ipad][irow][1]) * (xc2[ipad][irow][1] - x1);
            else
               sw1 = 0.5 * (y2 - yc2[ipad][irow][1]) * (xc1[ipad][irow][1] - x1);
            sw2         = (0.5 * (yc1[ipad][irow][0] + yc2[ipad][irow][0]) - y1) * wPad;
            SPads[ipad] = wPad * hPad - sw1 - sw2;
            continue;
         }
         // case g)
         if (indexPad[ipad][irow][0] == 1100 && indexPad[ipad][irow][1] == 1100) {
            if (xc1[ipad][irow][0] == x1) {
               sw1 = 0.5 * (y2 - yc1[ipad][irow][0]) * (xc2[ipad][irow][0] - x1);
               if (xc1[ipad][irow][1] == x1)
                  sw2 = 0.5 * (y2 - yc1[ipad][irow][1]) * (xc2[ipad][irow][1] - x1);
               else
                  sw2 = 0.5 * (y2 - yc2[ipad][irow][1]) * (xc1[ipad][irow][1] - x1);
            } else {
               sw1 = 0.5 * (y2 - yc2[ipad][irow][0]) * (xc1[ipad][irow][0] - x1);
               if (xc1[ipad][irow][1] == x1)
                  sw2 = 0.5 * (y2 - yc1[ipad][irow][1]) * (xc2[ipad][irow][1] - x1);
               else
                  sw2 = 0.5 * (y2 - yc2[ipad][irow][1]) * (xc1[ipad][irow][1] - x1);
            }
            SPads[ipad] = TMath::Abs(sw1 - sw2);
            continue;
         }
         // *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
         if (indexPad[ipad][irow][0] == 101 && indexPad[ipad][irow][1] == 101) {
            // case h)
            if (yc1[ipad][irow][0] == y1)
               if (yc1[ipad][irow][1] == y1)
                  sw1 = TMath::Abs(xc1[ipad][irow][0] - xc1[ipad][irow][1]);
               else
                  sw1 = TMath::Abs(xc1[ipad][irow][0] - xc2[ipad][irow][1]);
            else if (yc1[ipad][irow][1] == y1)
               sw1 = TMath::Abs(xc2[ipad][irow][0] - xc1[ipad][irow][1]);
            else
               sw1 = TMath::Abs(xc2[ipad][irow][0] - xc2[ipad][irow][1]);
            SPads[ipad] = sw1 * hPad;
            continue;
         }

         //                     base 0110
         // case b) or d)
         if (indexPad[ipad][irow][0] == 110 && indexPad[ipad][irow][1] == 0) {
            if (yc1[ipad][irow][0] > yc2[ipad][irow][0])
               sw1 = 0.5 * (y2 - yc2[ipad][irow][0]) * (x2 - xc1[ipad][irow][0]);
            else
               sw1 = 0.5 * (y2 - yc1[ipad][irow][0]) * (x2 - xc2[ipad][irow][0]);
            if (ipad < ip2 - 1)
               if (indexPad[ipad + 1][irow][1] != 0)
                  SPads[ipad] = sw1; // d)
               else
                  SPads[ipad] = wPad * hPad - sw1; // b)
            else if (indexPad[ipad - 1][irow][1] != 0)
               SPads[ipad] = wPad * hPad - sw1; // b)
            else
               SPads[ipad] = sw1; // d)
            continue;
         }
         if (indexPad[ipad][irow][1] == 110 && indexPad[ipad][irow][0] == 0) {
            if (yc1[ipad][irow][1] > yc2[ipad][irow][1])
               sw1 = 0.5 * (y2 - yc2[ipad][irow][1]) * (x2 - xc1[ipad][irow][1]);
            else
               sw1 = 0.5 * (y2 - yc1[ipad][irow][1]) * (x2 - xc2[ipad][irow][1]);
            if (ipad < ip2 - 1)
               if (indexPad[ipad + 1][irow][0] != 0)
                  SPads[ipad] = sw1; // d)
               else
                  SPads[ipad] = wPad * hPad - sw1; // b)
            else if (indexPad[ipad - 1][irow][0] != 0)
               SPads[ipad] = wPad * hPad - sw1; // b)
            else
               SPads[ipad] = sw1; // d)
            continue;
         }
         // case c)
         if (indexPad[ipad][irow][0] == 1010 && indexPad[ipad][irow][1] == 0) {
            if (irow < nRows - 1) {
               if (indexPad[ipad][irow + 1][1] != 0)
                  SPads[ipad] = (y2 - 0.5 * (yc1[ipad][irow][0] + yc2[ipad][irow][0])) * wPad;
               else
                  SPads[ipad] = (0.5 * (yc1[ipad][irow][0] + yc2[ipad][irow][0]) - y1) * hPad;
               continue;
            } else {
               if (indexPad[ipad - 1][irow][1] != 0)
                  SPads[ipad] = (0.5 * (yc1[ipad][irow][0] + yc2[ipad][irow][0]) - y1) * hPad;
               else
                  SPads[ipad] = (y2 - 0.5 * (yc1[ipad][irow][0] + yc2[ipad][irow][0])) * wPad;
               continue;
            }
         }
         if (indexPad[ipad][irow][1] == 1010 && indexPad[ipad][irow][0] == 0) {
            if (irow < nRows - 1) {
               if (indexPad[ipad][irow + 1][1] != 0)
                  SPads[ipad] = (y2 - 0.5 * (yc1[ipad][irow][1] + yc2[ipad][irow][1])) * wPad;
               else
                  SPads[ipad] = (0.5 * (yc1[ipad][irow][1] + yc2[ipad][irow][1]) - y1) * hPad;
               continue;
            } else {
               if (indexPad[ipad - 1][irow][0] != 0)
                  SPads[ipad] = (0.5 * (yc1[ipad][irow][1] + yc2[ipad][irow][1]) - y1) * hPad;
               else
                  SPads[ipad] = (y2 - 0.5 * (yc1[ipad][irow][1] + yc2[ipad][irow][1])) * wPad;
               continue;
            }
         }
         // case e)
         if (indexPad[ipad][irow][0] == 110 && indexPad[ipad][irow][1] == 1001) {
            if (yc1[ipad][irow][1] == y1) {
               sw1 = 0.5 * (yc2[ipad][irow][1] - y1) * (xc1[ipad][irow][1] - x1);
               if (yc1[ipad][irow][0] > yc2[ipad][irow][0])
                  sw2 = 0.5 * (x2 - xc1[ipad][irow][0]) * (y2 - yc2[ipad][irow][0]);
               else
                  sw2 = 0.5 * (x2 - xc2[ipad][irow][0]) * (y2 - yc1[ipad][irow][0]);
            } else {
               sw1 = 0.5 * (yc1[ipad][irow][1] - y1) * (xc2[ipad][irow][1] - x1);
               if (yc1[ipad][irow][0] > yc2[ipad][irow][0])
                  sw2 = 0.5 * (x2 - xc1[ipad][irow][0]) * (y2 - yc2[ipad][irow][0]);
               else
                  sw2 = 0.5 * (x2 - xc2[ipad][irow][0]) * (y2 - yc1[ipad][irow][0]);
            }
            SPads[ipad] = wPad * hPad - sw1 - sw2;
            continue;
         }
         if (indexPad[ipad][irow][1] == 110 && indexPad[ipad][irow][0] == 1001) {
            if (yc1[ipad][irow][0] == y1) {
               sw1 = 0.5 * (yc2[ipad][irow][0] - y1) * (xc1[ipad][irow][0] - x1);
               if (yc1[ipad][irow][1] > yc2[ipad][irow][1])
                  sw2 = 0.5 * (x2 - xc1[ipad][irow][1]) * (y2 - yc2[ipad][irow][1]);
               else
                  sw2 = 0.5 * (x2 - xc2[ipad][irow][1]) * (y2 - yc1[ipad][irow][1]);
            } else {
               sw1 = 0.5 * (yc1[ipad][irow][0] - y1) * (xc2[ipad][irow][0] - x1);
               if (yc1[ipad][irow][1] > yc2[ipad][irow][1])
                  sw2 = 0.5 * (x2 - xc1[ipad][irow][1]) * (y2 - yc2[ipad][irow][1]);
               else
                  sw2 = 0.5 * (x2 - xc2[ipad][irow][1]) * (y2 - yc1[ipad][irow][1]);
            }
            SPads[ipad] = wPad * hPad - sw1 - sw2;
            continue;
         }
         // case f)
         if (indexPad[ipad][irow][0] == 110 && indexPad[ipad][irow][1] == 101) {
            if (yc1[ipad][irow][0] > yc2[ipad][irow][0])
               sw1 = 0.5 * (x2 - xc1[ipad][irow][0]) * (y2 - yc2[ipad][irow][0]);
            else
               sw1 = 0.5 * (x2 - xc2[ipad][irow][0]) * (y2 - yc1[ipad][irow][0]);
            sw2         = (0.5 * (xc1[ipad][irow][1] + xc2[ipad][irow][1]) - x1) * hPad;
            SPads[ipad] = wPad * hPad - sw1 - sw2;
            continue;
         }
         if (indexPad[ipad][irow][1] == 110 && indexPad[ipad][irow][0] == 101) {
            if (yc1[ipad][irow][1] > yc2[ipad][irow][1])
               sw1 = 0.5 * (x2 - xc1[ipad][irow][1]) * (y2 - yc2[ipad][irow][1]);
            else
               sw1 = 0.5 * (x2 - xc2[ipad][irow][1]) * (y2 - yc1[ipad][irow][1]);
            sw2         = (0.5 * (xc1[ipad][irow][0] + xc2[ipad][irow][0]) - x1) * hPad;
            SPads[ipad] = wPad * hPad - sw1 - sw2;
            continue;
         }
         if (indexPad[ipad][irow][0] == 110 && indexPad[ipad][irow][1] == 1010) {
            if (yc1[ipad][irow][0] > yc2[ipad][irow][0])
               sw1 = 0.5 * (x2 - xc1[ipad][irow][0]) * (y2 - yc2[ipad][irow][0]);
            else
               sw1 = 0.5 * (x2 - xc2[ipad][irow][0]) * (y2 - yc1[ipad][irow][0]);
            sw2         = (0.5 * (yc1[ipad][irow][1] + yc2[ipad][irow][1]) - y1) * wPad;
            SPads[ipad] = wPad * hPad - sw1 - sw2;
            continue;
         }
         if (indexPad[ipad][irow][1] == 110 && indexPad[ipad][irow][0] == 1010) {
            if (yc1[ipad][irow][1] > yc2[ipad][irow][1])
               sw1 = 0.5 * (x2 - xc1[ipad][irow][1]) * (y2 - yc2[ipad][irow][1]);
            else
               sw1 = 0.5 * (x2 - xc2[ipad][irow][1]) * (y2 - yc1[ipad][irow][1]);
            sw2         = (0.5 * (yc1[ipad][irow][0] + yc2[ipad][irow][0]) - y1) * wPad;
            SPads[ipad] = wPad * hPad - sw1 - sw2;
            continue;
         }
         // case g)
         if (indexPad[ipad][irow][0] == 110 && indexPad[ipad][irow][1] == 110) {
            if (yc1[ipad][irow][0] > yc2[ipad][irow][0]) {
               sw1 = 0.5 * (x2 - xc2[ipad][irow][0]) * (y2 - yc1[ipad][irow][0]);
               if (yc1[ipad][irow][1] > yc2[ipad][irow][1])
                  sw2 = 0.5 * (x2 - xc1[ipad][irow][1]) * (y2 - yc2[ipad][irow][1]);
               else
                  sw2 = 0.5 * (x2 - xc2[ipad][irow][1]) * (y2 - yc1[ipad][irow][1]);
            } else {
               sw1 = 0.5 * (x2 - xc1[ipad][irow][0]) * (y2 - yc2[ipad][irow][0]);
               if (yc1[ipad][irow][1] > yc2[ipad][irow][1])
                  sw2 = 0.5 * (x2 - xc1[ipad][irow][1]) * (y2 - yc2[ipad][irow][1]);
               else
                  sw2 = 0.5 * (x2 - xc2[ipad][irow][1]) * (y2 - yc1[ipad][irow][1]);
            }
            SPads[ipad] = TMath::Abs(sw1 - sw2);
            continue;
         }
         // *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
         if (indexPad[ipad][irow][0] == 1010 && indexPad[ipad][irow][1] == 1010) {
            // case h)
            if (xc1[ipad][irow][0] == x1)
               if (xc1[ipad][irow][1] == x1)
                  sw1 = TMath::Abs(yc1[ipad][irow][0] - yc1[ipad][irow][1]);
               else
                  sw1 = TMath::Abs(yc1[ipad][irow][0] - yc2[ipad][irow][1]);
            else if (xc1[ipad][irow][1] == x1)
               sw1 = TMath::Abs(yc2[ipad][irow][0] - yc1[ipad][irow][1]);
            else
               sw1 = TMath::Abs(yc2[ipad][irow][0] - yc2[ipad][irow][1]);
            SPads[ipad] = sw1 * wPad;
            continue;
         }

         //                     base 0011
         // case b) or d)
         if (indexPad[ipad][irow][0] == 11 && indexPad[ipad][irow][1] == 0) {
            if (yc1[ipad][irow][0] < yc2[ipad][irow][0])
               sw1 = 0.5 * (x2 - xc1[ipad][irow][0]) * (yc2[ipad][irow][0] - y1);
            else
               sw1 = 0.5 * (x2 - xc2[ipad][irow][0]) * (yc1[ipad][irow][0] - y1);
            if (ipad < ip2 - 1)
               if (indexPad[ipad + 1][irow][1] != 0)
                  SPads[ipad] = sw1; // d)
               else
                  SPads[ipad] = wPad * hPad - sw1; // b)
            else if (indexPad[ipad - 1][irow][1] != 0)
               SPads[ipad] = wPad * hPad - sw1; // b)
            else
               SPads[ipad] = sw1; // d)
            continue;
         }
         if (indexPad[ipad][irow][1] == 11 && indexPad[ipad][irow][0] == 0) {
            if (yc1[ipad][irow][1] < yc2[ipad][irow][1])
               sw1 = 0.5 * (x2 - xc1[ipad][irow][1]) * (yc2[ipad][irow][1] - y1);
            else
               sw1 = 0.5 * (x2 - xc2[ipad][irow][1]) * (yc1[ipad][irow][1] - y1);
            if (ipad < ip2 - 1)
               if (indexPad[ipad + 1][irow][0] != 0)
                  SPads[ipad] = sw1; // d)
               else
                  SPads[ipad] = wPad * hPad - sw1; // b)
            else if (indexPad[ipad - 1][irow][0] != 0)
               SPads[ipad] = wPad * hPad - sw1; // b)
            else
               SPads[ipad] = sw1; // d)
            continue;
         }
         // case f)
         if (indexPad[ipad][irow][0] == 11 && indexPad[ipad][irow][1] == 101) {
            if (yc1[ipad][irow][0] == y1)
               sw1 = 0.5 * (x2 - xc1[ipad][irow][0]) * (yc2[ipad][irow][0] - y1);
            else
               sw1 = 0.5 * (x2 - xc2[ipad][irow][0]) * (yc1[ipad][irow][0] - y1);
            sw2         = (0.5 * (xc1[ipad][irow][1] + xc2[ipad][irow][1]) - x1) * hPad;
            SPads[ipad] = wPad * hPad - sw1 - sw2;
            continue;
         }
         if (indexPad[ipad][irow][1] == 0011 && indexPad[ipad][irow][0] == 0101) {
            if (yc1[ipad][irow][1] == y1)
               sw1 = 0.5 * (x2 - xc1[ipad][irow][1]) * (yc2[ipad][irow][1] - y1);
            else
               sw1 = 0.5 * (x2 - xc2[ipad][irow][1]) * (yc1[ipad][irow][1] - y1);
            sw2         = (0.5 * (xc1[ipad][irow][0] + xc2[ipad][irow][0]) - x1) * hPad;
            SPads[ipad] = wPad * hPad - sw1 - sw2;
            continue;
         }
         if (indexPad[ipad][irow][0] == 11 && indexPad[ipad][irow][1] == 1010) {
            if (yc1[ipad][irow][0] == y1)
               sw1 = 0.5 * (x2 - xc1[ipad][irow][0]) * (yc2[ipad][irow][0] - y1);
            else
               sw1 = 0.5 * (x2 - xc2[ipad][irow][0]) * (yc1[ipad][irow][0] - y1);
            sw2         = (y2 - 0.5 * (yc1[ipad][irow][1] + yc2[ipad][irow][1])) * wPad;
            SPads[ipad] = wPad * hPad - sw1 - sw2;
            continue;
         }
         if (indexPad[ipad][irow][1] == 11 && indexPad[ipad][irow][0] == 1010) {
            if (yc1[ipad][irow][1] == y1)
               sw1 = 0.5 * (x2 - xc1[ipad][irow][1]) * (yc2[ipad][irow][1] - y1);
            else
               sw1 = 0.5 * (x2 - xc2[ipad][irow][1]) * (yc1[ipad][irow][1] - y1);
            sw2         = (y2 - 0.5 * (yc1[ipad][irow][0] + yc2[ipad][irow][0])) * wPad;
            SPads[ipad] = wPad * hPad - sw1 - sw2;
            continue;
         }
         // case g)
         if (indexPad[ipad][irow][0] == 11 && indexPad[ipad][irow][1] == 11) {
            if (yc1[ipad][irow][0] == y1) {
               sw1 = 0.5 * (x2 - xc1[ipad][irow][0]) * (yc2[ipad][irow][0] - y1);
               if (yc1[ipad][irow][1] == y1)
                  sw2 = 0.5 * (x2 - xc1[ipad][irow][1]) * (yc2[ipad][irow][1] - y1);
               else
                  sw2 = 0.5 * (x2 - xc2[ipad][irow][1]) * (yc1[ipad][irow][1] - y1);
            } else {
               sw1 = 0.5 * (x2 - xc2[ipad][irow][0]) * (yc1[ipad][irow][0] - y1);
               if (yc1[ipad][irow][1] == y1)
                  sw2 = 0.5 * (x2 - xc1[ipad][irow][1]) * (yc2[ipad][irow][1] - y1);
               else
                  sw2 = 0.5 * (x2 - xc2[ipad][irow][1]) * (yc1[ipad][irow][1] - y1);
            }
            SPads[ipad] = TMath::Abs(sw1 - sw2);
            continue;
         }

         //                     base 1001
         // case b) or d)
         if (indexPad[ipad][irow][0] == 1001 && indexPad[ipad][irow][1] == 0) {
            if (xc1[ipad][irow][0] = x1)
               sw1 = 0.5 * (xc2[ipad][irow][0] - x1) * (yc1[ipad][irow][0] - y1);
            else
               sw1 = 0.5 * (xc1[ipad][irow][0] - x1) * (yc2[ipad][irow][0] - y1);
            if (ipad < ip2 - 1)
               if (indexPad[ipad + 1][irow][1] != 0)
                  SPads[ipad] = wPad * hPad - sw1; // b)
               else
                  SPads[ipad] = sw1; // d)
            else if (indexPad[ipad - 1][irow][1] != 0)
               SPads[ipad] = sw1; // d)
            else
               SPads[ipad] = wPad * hPad - sw1; // b)
            continue;
         }
         if (indexPad[ipad][irow][1] == 1001 && indexPad[ipad][irow][0] == 0) {
            if (xc1[ipad][irow][1] = x1)
               sw1 = 0.5 * (xc2[ipad][irow][1] - x1) * (yc1[ipad][irow][1] - y1);
            else
               sw1 = 0.5 * (xc1[ipad][irow][1] - x1) * (yc2[ipad][irow][1] - y1);
            if (ipad < ip2 - 1)
               if (indexPad[ipad + 1][irow][0] != 0)
                  SPads[ipad] = wPad * hPad - sw1; // b)
               else
                  SPads[ipad] = sw1; // d)
            else if (indexPad[ipad - 1][irow][0] != 0)
               SPads[ipad] = sw1; // d)
            else
               SPads[ipad] = wPad * hPad - sw1; // b)
            continue;
         }
         // case f)
         if (indexPad[ipad][irow][0] == 1001 && indexPad[ipad][irow][1] == 101) {
            if (xc1[ipad][irow][0] == x1)
               sw1 = 0.5 * (xc2[ipad][irow][0] - x1) * (yc1[ipad][irow][0] - y2);
            else
               sw1 = 0.5 * (xc1[ipad][irow][0] - x1) * (yc2[ipad][irow][0] - y2);
            sw2         = (x2 - 0.5 * (xc1[ipad][irow][1] + xc2[ipad][irow][1])) * hPad;
            SPads[ipad] = wPad * hPad - sw1 - sw2;
            continue;
         }
         if (indexPad[ipad][irow][1] == 1001 && indexPad[ipad][irow][0] == 101) {
            if (xc1[ipad][irow][1] == x1)
               sw1 = 0.5 * (xc2[ipad][irow][1] - x1) * (yc1[ipad][irow][1] - y2);
            else
               sw1 = 0.5 * (xc1[ipad][irow][1] - x1) * (yc2[ipad][irow][1] - y2);
            sw2         = (x2 - 0.5 * (xc1[ipad][irow][0] + xc2[ipad][irow][0])) * hPad;
            SPads[ipad] = wPad * hPad - sw1 - sw2;
            continue;
         }
         if (indexPad[ipad][irow][0] == 1001 && indexPad[ipad][irow][1] == 1010) {
            if (xc1[ipad][irow][0] == x1)
               sw1 = 0.5 * (xc2[ipad][irow][0] - x1) * (yc1[ipad][irow][0] - y2);
            else
               sw1 = 0.5 * (xc1[ipad][irow][0] - x1) * (yc2[ipad][irow][0] - y2);
            sw2         = (y2 - 0.5 * (yc1[ipad][irow][1] + yc2[ipad][irow][1])) * wPad;
            SPads[ipad] = wPad * hPad - sw1 - sw2;
            continue;
         }
         if (indexPad[ipad][irow][1] == 1001 && indexPad[ipad][irow][0] == 1010) {
            if (xc1[ipad][irow][1] == x1)
               sw1 = 0.5 * (xc2[ipad][irow][1] - x1) * (yc1[ipad][irow][1] - y2);
            else
               sw1 = 0.5 * (xc1[ipad][irow][1] - x1) * (yc2[ipad][irow][1] - y2);
            sw2         = (y2 - 0.5 * (yc1[ipad][irow][0] + yc2[ipad][irow][0])) * wPad;
            SPads[ipad] = wPad * hPad - sw1 - sw2;
            continue;
         }
         // case g)
         if (indexPad[ipad][irow][0] == 1001 && indexPad[ipad][irow][1] == 1001) {
            if (xc1[ipad][irow][0] == x1) {
               sw1 = 0.5 * (xc2[ipad][irow][0] - x1) * (yc1[ipad][irow][0] - y2);
               if (xc1[ipad][irow][1] == x1)
                  sw1 = 0.5 * (xc2[ipad][irow][1] - x1) * (yc1[ipad][irow][1] - y2);
               else
                  sw1 = 0.5 * (xc1[ipad][irow][1] - x1) * (yc2[ipad][irow][1] - y2);
            } else {
               sw1 = 0.5 * (xc1[ipad][irow][0] - x1) * (yc2[ipad][irow][0] - y2);
               if (xc1[ipad][irow][1] == x1)
                  sw1 = 0.5 * (xc2[ipad][irow][1] - x1) * (yc1[ipad][irow][1] - y2);
               else
                  sw1 = 0.5 * (xc1[ipad][irow][1] - x1) * (yc2[ipad][irow][1] - y2);
            }
            SPads[ipad] = TMath::Abs(sw1 - sw2);
            continue;
         }
         // *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
      } // end of pads cycle

      // for(Int_t p=0;p<ip2;p++)
      //  if(SPads[p]!=0) printf("S[%d,%d]=%f \n",p,irow,SPads[p]);
      //    ***************************************************
      //    ********       check pads in each row      ********
      //    ***************************************************
      Int_t iRoc = fTpcSecGeo->NofRoc(irow);
      Int_t iP2  = ip2 - 1;
      nclust += ClustersInRow(irow, ip1, iP2, iRoc, SPads);
      // printf("nclust=%d    irow,ip1,ip2= %d   %d-%d\n",nclust,irow,ip1,ip2);

   } // end of row  cycle
   // clear pad indexes for operationd with next sector
   for (Int_t i = 0; i < mxp; i++)
      for (Int_t j = 0; j < mxr; j++)
         for (Int_t k = 0; k < 3; k++) indexPad[i][j][k] = 0;
}

//______________________________________________________________________
// Get get helix point coordinates and sector number
Int_t TpcMiniMC::coordinates(Double_t t, Double_t R)
{
   Double_t phi;
   Int_t    sec;
   rBt = TVector3(xcB + R * TMath::Cos(phi_0B + signRotB * t), // BCS
                  ycB + R * TMath::Sin(phi_0B + signRotB * t), r0BeamB.Z() + helix_Pz * t);
   // printf("rBt(%f,%f,%f) rcB(%f,%f) phi_0B=%f\n",rBt.X(),rBt.Y(),rBt.Z(),
   //         xcB,ycB,phi_0B*r2d);
   rGt = B2G * rBt;
   if (TMath::Abs(rGt.X()) > epsMin || TMath::Abs(rGt.Y()) > epsMin) { // GCS
      if (rGt.Z() > 0) {
         phi = TMath::ATan2(rGt.Y(), rGt.X());
         if (TMath::Abs(phi) < phi_shift)
            sec = 0;
         else {
            if (phi < 0) phi += twopi;
            sec = Int_t((phi + phi_shift) / phi_sec);
         }
      } else {
         phi = TMath::ATan2(rGt.Y(), -rGt.X());
         if (TMath::Abs(phi) < phi_shift)
            sec = 0;
         else {
            if (phi < 0) phi += twopi;
            sec = Int_t((phi + phi_shift) / phi_sec) + 12;
         }
      }
   } else {
      sec = 0;
   }
   if (sec != iSector) {
      fTpcSecGeo->SectorTransformation(sec, R0shift, loc2glo);
      e_3 = TVector3(loc2glo.XZ(), loc2glo.YZ(), loc2glo.ZZ());
      fTpcSecGeo->SectorBackTransformation(sec, G0shift, glo2loc);
   }
   // printf("rBt(%f,%f,%f) rcB(%f,%f) phi_0B=%f  sec=%d\n",rBt.X(),rBt.Y(),rBt.Z(),
   //         xcB,ycB,phi_0B*r2d,sec);
   Double_t t_s   = (e_3 * (R0shift - rGt)) / (e_3 * eEf); // t to sector plane point
   TVector3 rGt_s = rGt + t_s * eEf;                       // helix point on the sector plane (GCS)
   rLt            = G0shift + glo2loc * rGt_s;
   return sec;
}

//______________________________________________________________________

// simulate z pozition
Double_t TpcMiniMC::simHelixZpos(TVector2 LCxy)
{

   // a) lower the perpendicular from LCxy to the track center line
   //    to the track of the particle.
   // This will be the average value of the Z coordinate measurement
   TVector3 SecL = TVector3(LCxy.X(), LCxy.Y(), 0.);
   TVector3 SecG = R0shift + loc2glo * SecL;
   TVector3 SecB = G2B * SecG;
   // printf("simHelixZpos: rL(%f,%f)  SecG(%f,%f,%f) SecB(%f,%f,%f)\n",
   // LCxy.X(),LCxy.Y(),SecG.X(),SecG.Y(),SecG.Z(),SecB.X(),SecB.Y(),SecB.Z());
   //  find the t helix parameter which corresponds to the cluster
   Double_t phi = TMath::ATan2(SecB.Y() - ycB, SecB.X() - xcB);
   Double_t t   = (q_sign * (phi_0B - phi) > 0) ? abs(phi_0B - phi) : twopi - abs(phi_0B - phi);
   // printf("simHelixZpos:  phi=%f  t=%f\n",
   // TMath::ATan2(SecB.Y()-ycB,SecB.X()-xcB)*r2d,t*r2d);
   //  find particle coordinates in BCS
   TVector3 particleB = TVector3(xcB + rH * TMath::Cos(phi_0B + signRotB * t), // BCS
                                 ycB + rH * TMath::Sin(phi_0B + signRotB * t), r0BeamB.Z() + helix_Pz * t);
   // find particle coordinates in GCS
   TVector3 particleG = B2G * particleB;
   // find the distance (D) from the particle to the sector plane
   Double_t t_s = (e_3 * (R0shift - particleG)) / (e_3 * eEf); // t to sector plane point
   TVector3 D   = t_s * eEf;                                   // distance to the track point from the cluster
   Double_t z   = fRandom->Gaus(D.Mag(), D.Mag() * sig_Z_pos);
   if (z < 0) {
      z = 0;
   } else {
      if (z > fL_svol) {
         z = fL_svol - sZ_min;
      }
   }
   // printf("simHelixZpos: particleB(%f,%f,%f) d=%f z=%f\n",
   // particleB.X(),particleB.Y(),particleB.Z(),D.Mag(),z);
   return z;
}

//______________________________________________________________________
// Chi2 for the clusters fit by the helix
Double_t TpcMiniMC_Chi2Helix(const Double_t *par)
{

   double chisq = 0;
   double delta;
   Int_t  nhits = MpdTpcMiniMC_Nhits;
   double phi0  = par[4];
   int    q     = MpdTpcMiniMC_qHelix;
   for (int i = 0; i < nhits; i++) {
      double X   = MpdTpcMiniMC_XBHits[i];
      double Y   = MpdTpcMiniMC_YBHits[i];
      double Z   = MpdTpcMiniMC_ZBHits[i];
      double W   = MpdTpcMiniMC_WHits[i];
      double phi = TMath::ATan2(Y - par[1], X - par[0]);
      double t;
      if (q > 0)
         t = (phi > phi0) ? TMath::TwoPi() - q * (phi - phi0) : -q * (phi - phi0);
      else
         t = (phi > phi0) ? -q * (phi - phi0) : TMath::TwoPi() - q * (phi - phi0);
      double dx = X - (par[0] + par[3] * TMath::Cos(par[4] - MpdTpcMiniMC_qHelix * t));
      double dy = Y - (par[1] + par[3] * TMath::Sin(par[4] - MpdTpcMiniMC_qHelix * t));
      double dz = Z - (par[2] + par[5] * t);
      delta     = TMath::Sqrt(dx * dx + dy * dy + dz * dz) / W;
      chisq += delta * delta;
   }
   return chisq;
}
// printf("   chi2[%d]=%f\n",i,chisq);

//______________________________________________________________________
// Reconstruction of the helix parameters
Int_t TpcMiniMC::fitHelix()
{
   /*
     in the BCS the helix equation:
       x=p[0]+p[3]*cos(p[4]+q*t)
       y=p[1]+p[3]*sin(p[4]+q*t)
       z=p[2]+p[5]*t
   */
   /*trackPout[0]=trackPin[0];
   trackPout[1]=trackPin[1];
   trackPout[2]=trackPin[2];
   trackPout[3]=trackPin[3];
   trackPout[4]=trackPin[4];
   trackPout[5]=trackPin[5];
   for(Int_t k=0;k<6;k++) trackDiff[k]=trackPout[k]-trackPin[k];
   return 0;*/
   // initial approach
   /*for(Int_t i=0;i<MpdTpcMiniMC_Nhits;i++) {
   printf("rG(%f,%f,%f) rB(%f,%f,%f)  s=%d\n",
   MpdTpcMiniMC_XHits[i],MpdTpcMiniMC_YHits[i],MpdTpcMiniMC_ZHits[i],
   MpdTpcMiniMC_XBHits[i],MpdTpcMiniMC_YBHits[i],MpdTpcMiniMC_ZBHits[i],
   MpdTpcMiniMC_iSec[i]);}*/
   Int_t    imidle = MpdTpcMiniMC_Nhits / 2;
   Double_t x1, y1, z1, a1, b1, c1, x2, y2, z2, a2, b2, c2, x3, y3, z3, xc, yc, xc1223, yc1223, xc1213, yc1213, xc1323,
      yc1323, d, t, t1, t2, t3, q;
   x1 = MpdTpcMiniMC_XBHits[0];
   y1 = MpdTpcMiniMC_YBHits[0];
   z1 = MpdTpcMiniMC_ZBHits[0];
   x2 = MpdTpcMiniMC_XBHits[imidle];
   y2 = MpdTpcMiniMC_YBHits[imidle];
   x3 = MpdTpcMiniMC_XBHits[MpdTpcMiniMC_Nhits - 1];
   y3 = MpdTpcMiniMC_YBHits[MpdTpcMiniMC_Nhits - 1];
   z3 = MpdTpcMiniMC_ZBHits[MpdTpcMiniMC_Nhits - 1];

   a1 = x2 - x1;
   b1 = y2 - y1;
   a2 = x3 - x1;
   b2 = y3 - y1;
   c1 = 0.5 * (x2 * x2 - x1 * x1 + y2 * y2 - y1 * y1);
   c2 = 0.5 * (x3 * x3 - x1 * x1 + y3 * y3 - y1 * y1);
   d  = a1 * b2 - a2 * b1;
   if (TMath::Abs(d) > epsMin) { // crossing of 2 perpendiculars
      xc = (c1 * b2 - c2 * b1) / d;
      yc = (a1 * c2 - a2 * c1) / d;
   } else { // 3 points lie on a line
      TVector3 v  = TVector3(x3 - x1, y3 - y1, 0);
      TVector3 B  = TVector3(0, 0, 1);
      TVector3 rc = v.Cross(B);
      xc          = 0.5 * (x3 + x1) + q * 100000 * rc.X() / rc.Mag();
      yc          = 0.5 * (y3 + y1) + q * 100000 * rc.Y() / rc.Mag();
   }
   //    fitCenter(xc,yc,xc1213,yc1213);
   //      xc=xc1213;
   //      yc=yc1213;

   // determ the particle sign
   q  = 0;
   t1 = TMath::ATan2(y1 - yc, x1 - xc);
   for (Int_t i = 1; i < MpdTpcMiniMC_Nhits; i++) {
      t2 = TMath::ATan2(MpdTpcMiniMC_YBHits[i] - yc, MpdTpcMiniMC_XBHits[i] - xc);
      t  = t2 - t1;
      if (TMath::Abs(t) > pi)
         if (t1 < 0)
            q += twopi - (t2 - t1);
         else
            q += (t1 - t2) - twopi;
      else
         q += t;
      t1 = t2;
   }
   MpdTpcMiniMC_qHelix = (q > 0) ? -1 : 1; // sign of the particle
   // printf("r1(%f,%f,%f) r2(%f,%f,%f) r3(%f,%f,%f)\n",x1,y1,z1,x2,y2,z2,x3,y3,z3);
   //  Choose method upon creation between:
   //  kMigrad, kSimplex, kCombined,
   //  kScan, kFumili
   static ROOT::Minuit2::Minuit2Minimizer minFit(ROOT::Minuit2::kMigrad);

   ROOT::Math::Functor fcn(&TpcMiniMC_Chi2Helix, 6);
   double              step[6] = {1, 1, 1, 1, 0.1, 1};
   double              limU[6] = {100000, 100000, 100000, 100000, TMath::Pi(), 100000};
   double              limL[6] = {-100000, -100000, -100000, 0, -TMath::Pi(), -100000};
   double              par[6];

   par[0] = xc;                             //
   par[1] = yc;                             //  - coordinates of the muon start point
   par[2] = 0;                              //
   par[3] = TMath::Sqrt(xc * xc + yc * yc); // radius of the circle in XY plane of BCS
   par[4] = TMath::ATan2(-yc, -xc);         // initial muon angle (XY plane of BCS)
   par[5] = (z3 - z1) / TMath::Abs(q);

   // printf("fitHelix: par(%8.1f %8.1f %4.1f %8.1f %6.1f %8.1f)  q=%d\n",
   // par[0],par[1],par[2],par[3],r2d*par[4],par[5],MpdTpcMiniMC_qHelix);

   minFit.SetFunction(fcn);

   // Set the free variables to be minimized!
   minFit.SetLimitedVariable(0, "x0", par[0], step[0], limL[0], limU[0]);
   minFit.SetLimitedVariable(1, "y0", par[1], step[1], limL[1], limU[1]);
   minFit.SetLimitedVariable(2, "z0", par[2], step[2], limL[2], limU[2]);
   minFit.SetLimitedVariable(3, "R", par[3], step[3], limL[3], limU[3]);
   minFit.SetLimitedVariable(4, "phi", par[4], step[4], limL[4], limU[4]);
   minFit.SetLimitedVariable(5, "Pz", par[5], step[5], limL[5], limU[5]);

   Bool_t res = minFit.Minimize();
   //    if(res) {
   if (minFit.MinValue() / MpdTpcMiniMC_Nhits > EvalpHit) return 1;
   //      minFit.PrintResults();
   const double *xs = minFit.X();
   trackPout[0]     = xs[0];
   trackPout[1]     = xs[1];
   trackPout[2]     = xs[2];
   trackPout[3]     = MpdTpcMiniMC_qHelix * xs[3];
   trackPout[4]     = xs[4];
   trackPout[5]     = xs[5];
   for (Int_t k = 0; k < 6; k++) trackDiff[k] = trackPout[k] - trackPin[k];
   return 0;
   //    }
   //    else return 1;
}

//______________________________________________________________________
// put chi2 values
void TpcMiniMC::putChi2H()
{
   Double_t chisq = 0;
   Double_t delta;
   Int_t    nhits = MpdTpcMiniMC_Nhits;
   Double_t phi0  = trackPout[4];
   Double_t R     = TMath::Abs(trackPout[3]);
   Int_t    q     = trackPout[3] / R;
   /*printf("putChi2: trackPin=(%f,%f,%f,%f,%f,%f,)\n",
   trackPin[0],trackPin[1],trackPin[2],trackPin[3],trackPin[4],trackPin[5]);
   printf("putChi2: trackPout=(%f,%f,%f,%f,%f,%f,)\n",
   trackPout[0],trackPout[1],trackPout[2],trackPout[3],trackPout[4],trackPout[5]);*/
   for (int i = 0; i < nhits; i++) {
      Double_t X   = MpdTpcMiniMC_XBHits[i];
      Double_t Y   = MpdTpcMiniMC_YBHits[i];
      Double_t Z   = MpdTpcMiniMC_ZBHits[i];
      Double_t W   = MpdTpcMiniMC_WHits[i];
      Double_t phi = TMath::ATan2(Y - trackPout[1], X - trackPout[0]);
      Double_t t;
      if (q > 0)
         t = (phi > phi0) ? twopi - q * (phi - phi0) : -q * (phi - phi0);
      else
         t = (phi > phi0) ? -q * (phi - phi0) : twopi - q * (phi - phi0);
      Double_t dx = X - (trackPout[0] + R * TMath::Cos(phi0 - q * t));
      Double_t dy = Y - (trackPout[1] + R * TMath::Sin(phi0 - q * t));
      Double_t dz = Z - (trackPout[2] + trackPout[5] * t);
      delta       = TMath::Sqrt(dx * dx + dy * dy + dz * dz) / W;
      chisq += delta * delta;
      MpdTpcMiniMC_chi2s[i] = delta * delta;
      chisq += delta * delta;
      fNhitSec[MpdTpcMiniMC_iSec[i]]++;
      fSecChi2[MpdTpcMiniMC_iSec[i]] += delta * delta;
   }
   Chi2Fit = chisq / MpdTpcMiniMC_Nhits;
}
/*printf("putChi2: t=%f  secec[%d]=%d delta=%f  chi2=%f\n",t,i,
MpdTpcMiniMC_iSec[i],delta*delta,chisq);
printf("putChi2:   _r_i(%f,%f,%f)  r(%f,%f,%f \n",
X,Y,Z,trackPout[0]+R*TMath::Cos(phi0-q*t),trackPout[1]+R*TMath::Sin(phi0-q*t),
trackPout[2]+trackPout[5]*t);*/
// printf("putChi2: MpdTpcMiniMC_iSec[%d] delta^2=%f  chi2=%f\n",i,delta*delta,chisq);
//       if(Chi2Fit>100)
// printf("putChi2: hits=%d Chi2Fit per hit = %f\n",MpdTpcMiniMC_Nhits,Chi2Fit);

//______________________________________________________________________
// Get a ray of the TPC Laser System
void TpcMiniMC::GetLaserRay()
{
   const Double_t Fieldcage_shift_deg = 8.;
   const Double_t R = 123 - 3. / TMath::Cos(15. * TMath::DegToRad()); // Membrane_outer_holder_R_edge -
                                                                      // Fieldcage_Pin_R[1] /
                                                                      // TMath::Cos(Section_rad_step / 2.)
   const Double_t spread_ang = 11.;
   const Double_t pos_offset = 30.;
   Double_t       angle, x, y, phi, theta, z;

   // before first call ls1=true and ls2=false,
   if (ls1) { // check 1st ray section indexes
      if (a_ls1 > a12_ls) {
         a_ls1 = a11_ls;
         i_ls1++;
         if (i_ls1 >= i2_ls) {
            i_ls1 = i1_ls;
            j_ls1++;
            if (j_ls1 >= j2_ls) { // loop on 1st ray section is completed
               j_ls1 = j1_ls;
               ls1   = false;
               ls2   = true;
               a_ls2 = a21_ls;
               i_ls2 = i1_ls;
               j_ls2 = j1_ls;
            }
         }
      }
   }
   // printf("GetLaserRay:  1 1-%d(%d,%d,%d)   2-%d(%d,%d,%d)\n",
   // ls1,j_ls1,i_ls1,a_ls1,ls2,j_ls2,i_ls2,a_ls2);
   if (ls1) { // loop on 1st ray section is not completed
      angle = i_ls1 * 90. + Fieldcage_shift_deg;
      x     = R * TMath::Sin(angle * TMath::DegToRad());
      y     = R * TMath::Cos(angle * TMath::DegToRad());
      phi   = 360. - 90. * (i_ls1 + 1) - Fieldcage_shift_deg + spread_ang * a_ls1;
      theta = 90.;
      z     = j_ls1 < 0 ? j_ls1 * pos_offset : (j_ls1 + 1) * pos_offset;
      a_ls1++;
   }

   if (ls2) { // check 1st ray section indexes
      if (a_ls2 >= a22_ls) {
         a_ls2 = a21_ls;
         i_ls2++;
         if (i_ls2 >= i2_ls) {
            i_ls2 = i1_ls;
            j_ls2++;
            if (j_ls2 >= j2_ls) { // loop on 1st ray section is completed
               j_ls2 = j1_ls;
               ls2   = false;
               ls1   = true;
               a_ls1 = a11_ls;
               i_ls1 = i1_ls;
               j_ls1 = j1_ls;
            }
         }
      }
   }
   // printf("GetLaserRay:  2 1-%d(%d,%d,%d)   2-%d(%d,%d,%d)\n",
   // ls1,j_ls1,i_ls1,a_ls1,ls2,j_ls2,i_ls2,a_ls2);
   if (ls2) { // loop on 2nd ray section is not completed
      angle = i_ls2 * 90. + Fieldcage_shift_deg;
      x     = R * TMath::Sin(angle * TMath::DegToRad());
      y     = R * TMath::Cos(angle * TMath::DegToRad());
      phi   = 360. - 90. * (i_ls2 + 1) - Fieldcage_shift_deg +
            (a_ls2 < 0 ? (a_ls2 - 1) * spread_ang : (a_ls2 + 2) * spread_ang);
      theta = 90.;
      z     = j_ls2 < 0 ? j_ls2 * pos_offset : (j_ls2 + 1) * pos_offset;
      a_ls2++;
   }
   setMuon(x, y, z, theta * d2r, phi * d2r);
   // printf("GetLaserRay: R(%f,%f,%f) theta,phi(,%f,%f)\n",x,y,z,theta,phi);
}

//______________________________________________________________________
// Set single muon parameters
void TpcMiniMC::setMuon(Double_t x, Double_t y, Double_t z, Double_t theta, Double_t phi)
{
   r0Beam      = TVector3(x, y, z); // point inside TPC
   pLas        = TVector3(TMath::Sin(theta) * TMath::Cos(phi), TMath::Sin(theta) * TMath::Sin(phi),
                          TMath::Cos(theta)); // direction
   trackPin[0] = r0Beam.X();
   trackPin[1] = r0Beam.Y();
   trackPin[2] = r0Beam.Z();
   trackPin[3] = TMath::ACos(pLas.Z());
   trackPin[4] = TMath::ATan2(pLas.Y(), pLas.X());
}

//______________________________________________________________________
// Get random muon from the point
void TpcMiniMC::GetBeamMuon()
{
   Double_t theta1, theta2, theta, phi;
   // muon momentum
   theta1 = minMuonTheta;
   theta2 = TMath::Pi() - minMuonTheta;
   if (chambers == 1) {
      theta2 = TMath::Pi() / 2;
   } else {
      if (chambers == 2) {
         theta1 = TMath::Pi() / 2;
      } else {
         if (chambers != 0) {
            cout << "GetBeamMuon: Gno the mode chambers=" << chambers << endl;
            return;
         }
      }
   }
   theta = theta1 + (theta2 - theta1) * fRandom->Rndm();
   phi   = TMath::TwoPi() * fRandom->Rndm();
   pLas  = TVector3(TMath::Sin(theta) * TMath::Cos(phi), TMath::Sin(theta) * TMath::Sin(phi), TMath::Cos(theta));
}
//______________________________________________________________________
// Chi2 for the clusters fit by the helix
Double_t MpdTpcMiniMC_Chi2Center(const Double_t *par)
{

   double chisq = 0;
   double delta;
   Int_t  nhits = MpdTpcMiniMC_Nhits;
   for (int i = 0; i < nhits; i++) {
      double X  = MpdTpcMiniMC_XBHits[i];
      double Y  = MpdTpcMiniMC_YBHits[i];
      double dx = X - par[0];
      double dy = Y - par[1];
      chisq += dx * dx + dy * dy;
   }
   return chisq;
}

//______________________________________________________________________
// look for the helix center
Int_t TpcMiniMC::fitCenter(Double_t x0, Double_t y0, Double_t &xc, Double_t &yc)
{
   /*
     in the BCS the helix circle equation:
       (x-xc)**2+(y-yc)**2=R**2
   */
   // Choose method upon creation between:
   // kMigrad, kSimplex, kCombined,
   // kScan, kFumili
   static ROOT::Minuit2::Minuit2Minimizer minFit(ROOT::Minuit2::kMigrad);

   ROOT::Math::Functor fcn(&MpdTpcMiniMC_Chi2Center, 2);
   double              step[2] = {1, 1};
   double              limU[2] = {100000, 100000};
   double              limL[2] = {-100000, -100000};
   double              par[2];

   //    Int_t imidle=MpdTpcMiniMC_Nhits/2;
   // initial approach
   par[0] = x0;
   par[1] = y0;

   minFit.SetFunction(fcn);

   // Set the free variables to be minimized!
   minFit.SetLimitedVariable(0, "xc", par[0], step[0], limL[0], limU[0]);
   minFit.SetLimitedVariable(1, "yc", par[1], step[1], limL[1], limU[1]);

   minFit.Minimize();
   cout << "    FitMinResults\n";
   minFit.PrintResults();
   const double *xs = minFit.X();
   xc               = xs[0];
   yc               = xs[1];
   return 0;
}
