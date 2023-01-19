///
/// \brief MiniMC for TPC alignment tests
///
/// \author Valentin Kuzmin, SINP of Moscow State University

#include "TpcMissAlignment.h"
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

// Constructor difines all parameters
TpcMissAlignment::TpcMissAlignment(BaseTpcSectorGeo &tpcGeo, Double_t *B_mf, Bool_t Beam_point, Double_t *E_ef,
                                   Double_t *P_las, Double_t D_track, Double_t *Min_padS, Double_t *Max_padS,
                                   Double_t *Sig_padS, Double_t Sig_timeBins, Double_t Sig_Z_pos, Double_t *R0_beam,
                                   Double_t K_par, Double_t P_min, Double_t P_max, Int_t Chambers, Int_t MaxPcluster,
                                   Int_t MinHits, UInt_t Seed, TString *ARealIn, TString *AUsedIn, TString *OutDir,
                                   TString *OutFile = new TString(" "))
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
   // printf("MpdTpcMiniMC: k_par=%d\n",k_par);
   p_min       = P_min;
   p_max       = P_max;
   maxPcluster = MaxPcluster;
   d_track     = D_track;
   fSeed       = Seed;
   AReal       = ARealIn;
   AUsed       = AUsedIn;
   outDir      = OutDir;
   outFile     = OutFile;
   if (0 == 1) {
      printf("================== TcpMissAlignment parameters ====================\n");
      printf("B_mf(%f,%f,%f) BP=%d E_ef(%f,%f,%f) P_las(%f,%f,%f)\n", B_mf[0], B_mf[1], B_mf[2], Beam_point, E_ef[0],
             E_ef[1], E_ef[2], P_las[0], P_las[2], P_las[2]);
      printf("D_tr=%f Min_padS=(%f,%f) Max_pad=(%f,%f)\nSig_padS=(%f,%f) Sig_timeBins=%f Sig_Z_pos=%f\n", D_track,
             Min_padS[0], Min_padS[1], Max_padS[0], Max_padS[1], Sig_padS[0], Sig_padS[1], Sig_timeBins, Sig_Z_pos);
      printf("R0_beam(%f,%f,%f) K_par=%f P_min=%f P_max=%f\n", R0_beam[0], R0_beam[1], R0_beam[2], K_par, P_min, P_max);
      printf("Chambers=%d MaxPcluster=%d MinHits=%d Seed=%d\n", Chambers, MaxPcluster, MinHits, Seed);
      cout << "ARealIn=" << ARealIn->Data() << endl;
      cout << "AUsedIn=" << AUsedIn->Data() << endl;
      cout << "OutDir=" << OutDir->Data() << endl;
      cout << "OutFile=" << OutFile->Data() << endl;
      //  =%s\nAUsedIn=%s\nOutDir=%s\nOutFile=%s\n",
      //  ARealIn,AUsedIn,OutDir,OutFile);
      printf("===================================================================\n");
   }
   Init(tpcGeo);
}

// Constructor 3 difines reads parameters from the file
TpcMissAlignment::TpcMissAlignment(BaseTpcSectorGeo &tpcGeo, UInt_t Seed, TString *ParFile, TString *ARealIn,
                                   TString *AUsedIn, TString *OutDir, TString *OutFile)
{
   AReal   = ARealIn;
   AUsed   = AUsedIn;
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
void TpcMissAlignment::Init(BaseTpcSectorGeo &tpcGeo)
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
   if (!MFieldMode) {
      Double_t mod_p = sqrt(p_las[0] * p_las[0] + p_las[1] * p_las[1] + p_las[2] * p_las[2]);
      if (!MFieldMode)
         for (Int_t i = 0; i < 3; i++) p_las[i] = p_las[i] / mod_p;
   }
   // - - - - - - -
   fTpcSecGeo = dynamic_cast<TpcSectorGeoAlignmentVAK *>(&tpcGeo);
   if (!fTpcSecGeo)
      Fatal("TpcMissAlignment::Init", " !!! Wrong geometry type !!! "); // Zero alignment
                                                                        // - - - - - - -
   // read the real alignment
   const char *aInReal = AReal->Data();
   alignmentR          = new TFile(aInReal, "READ");
   ti_alignmentR       = (TTree *)alignmentR->Get("alignment");
   ti_alignmentR->SetBranchAddress("R0_A", &R0_AR);
   ti_alignmentR->SetBranchAddress("alpha_A", &alpha_AR);
   ti_alignmentR->SetBranchAddress("beta_A", &beta_AR);
   ti_alignmentR->SetBranchAddress("gamma_A", &gamma_AR);
   ti_alignmentR->GetEvent(0);
   delete ti_alignmentR;
   alignmentR->Close();
   fTpcSecGeo->TpcSectorGeoA(R0_AR, alpha_AR, beta_AR, gamma_AR);
   for (Int_t s = 0; s < 24; s++) {
      fTpcSecGeo->SectorBackTransformation(s, rG0_A[s], r_glo2loc[s]);
      fTpcSecGeo->SectorTransformation(s, rR0_A[s], r_loc2glo[s]);
   }
   // - - - - - - -
   // read the USED alignment
   const char *aInUsed = AUsed->Data();
   alignmentU          = new TFile(aInUsed, "READ");
   ti_alignmentU       = (TTree *)alignmentU->Get("alignment");
   ti_alignmentU->SetBranchAddress("R0_A", &R0_AU);
   ti_alignmentU->SetBranchAddress("alpha_A", &alpha_AU);
   ti_alignmentU->SetBranchAddress("beta_A", &beta_AU);
   ti_alignmentU->SetBranchAddress("gamma_A", &gamma_AU);
   ti_alignmentU->GetEvent(0);
   delete ti_alignmentU;
   alignmentU->Close();
   fTpcSecGeo->TpcSectorGeoA(R0_AU, alpha_AU, beta_AU, gamma_AU);
   for (Int_t s = 0; s < 24; s++) {
      fTpcSecGeo->SectorBackTransformation(s, uG0_A[s], u_glo2loc[s]);
      fTpcSecGeo->SectorTransformation(s, uR0_A[s], u_loc2glo[s]);
   }

   fL_svol = fTpcSecGeo->zPads().X();
   fZmin   = fTpcSecGeo->GetZmin();
   fZmax   = fTpcSecGeo->GetZmax();

   sZ_min = fZmax; // using in Helix z calculation

   bMf  = TVector3(b_mf[0], b_mf[1], b_mf[2]);
   bMag = bMf.Mag();

   eEf    = TVector3(e_ef[0], e_ef[1], e_ef[2]);
   r0Beam = TVector3(r0_beam[0], r0_beam[1], r0_beam[2]);
   pLas   = TVector3(p_las[0], p_las[1], p_las[2]);

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
   Double_t theta_B = TMath::ACos(bMf.Z() / bMag);
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
      Ssiz_t  ld     = AReal->Last('/');
      Ssiz_t  ln     = AReal->Last('.');
      TString fName0 = AReal->Data();
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
   t_parameters->Fill();
   t_parameters->Write();
   delete t_parameters;

   to_alignmentR = new TTree("real_alignment", "alignment_parameters");
   //  to_alignmentR = new TTree("alignment","alignment_parameters");
   to_alignmentR->Branch("R0_A", &R0_AR, "R0_AR[3][24]/D");
   to_alignmentR->Branch("alpha_A", &alpha_AR, "alpha_AR[24]/D");
   to_alignmentR->Branch("beta_A", &beta_AR, "beta_AR[24]/D");
   to_alignmentR->Branch("gamma_A", &gamma_AR, "gamma_AR[24]/D");
   to_alignmentR->Fill();
   to_alignmentR->Write();
   delete to_alignmentR;

   to_alignmentU = new TTree("alignment", "alignment_parameters");
   //  to_alignmentU  = new TTree("used_alignment","alignment_parameters");
   to_alignmentU->Branch("R0_A", &R0_AU, "R0_AU[3][24]/D");
   to_alignmentU->Branch("alpha_A", &alpha_AU, "alpha_AU[24]/D");
   to_alignmentU->Branch("beta_A", &beta_AU, "beta_AU[24]/D");
   to_alignmentU->Branch("gamma_A", &gamma_AU, "gamma_AU[24]/D");
   to_alignmentU->Fill();
   to_alignmentU->Write();
   delete to_alignmentU;

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
void TpcMissAlignment::simulate(Int_t nEvents)
{
   if (printL > 0) {
      cout << "simulate: simulation with Alignment "
           << "from the file: \n"
           << AReal->Data() << endl
           << "simulate: reconstruction with Alignment "
           << "from the file: \n"
           << AUsed->Data() << endl;
      cout << "simulate: " << nEvents << " with EvalCut=" << EvalpHit << endl;
   }
   fNevents       = nEvents;
   simEvt         = 0;
   Nreconstructed = 0;
   for (Int_t is = 0; is < 24; is++) {
      fNhitSec[is] = 0;
      fSecChi2[is] = 0;
   }
   while (Nreconstructed < fNevents) {
      simEvt++;
      Int_t attempts = miniEvent();
      // printf("simulate:  event=%d == attempts=%d nHits=%d =\n",simEvt,attempts,MpdTpcMiniMC_Nhits);
      if (attempts > 0 && MpdTpcMiniMC_Nhits >= fminHits) {
         if (printL > 1) {
            if (simEvt % printL == 0)
               cout << "MpdTpcMiniMC: event=" << simEvt << " simulated on the " << attempts
                    << " (th) attempt reconstructed=" << Nreconstructed << endl;
         }
         if (!MFieldMode) {
            if (fitLine() != 0) continue;
            // printf("   simulate: if(fitLine()==0)) =============\n");
            putChi2L();
         } else {
            // printf("   simulate: call fitHelix(%d) =============\n",simEvt);
            if (fitHelix() != 0) continue;
            putChi2H();
         }
         Nreconstructed++;
         // printf("   simulate: call putChi2(%d) =============\n",MpdTpcMiniMC_Nhits);
         t_data->Fill();
      } else {
         if (attempts < 0) {
            cout << "MpdTpcMiniMC: event=" << simEvt + 1 << " was not simulated on the " << attempts << " (th) attempt"
                 << endl;
            return;
         }
      }
      /*      if(simEvt>fNevents&&simEvt>3*Nreconstructed) {
              printf("STOP: very low efficienty of %6.2f at %d simulated events\n",
              float_t(Nreconstructed)/float_t(simEvt)*100,simEvt);
              return;
            }*/
   }
   t_data->Write();
   printf("******  %d events reconstucted of %d simulated (%%%4.1f) ******\n", Nreconstructed, simEvt,
          float_t(Nreconstructed) / float_t(simEvt) * 100);
   fSeed = fRandom->GetSeed();

   for (Int_t is = 0; is < 24; is++) {
      fNhits += fNhitSec[is];
      fSumChi2 += fSecChi2[is];
   }
   if (fNhits > 0) fSumChi2 /= fNhits;
   tout_fullChi2->Fill();
   tout_fullChi2->Write();
   delete tout_fullChi2;
   outmMCfile->Close();
   if (printL > 0) cout << " =============  FinalSeed = " << fSeed << endl;
}

//______________________________________________________________________
// an event generation
Int_t TpcMissAlignment::miniEvent()
{
   // printf("miniEvent: event=%d ================================\n",Nreconstructed);
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
                     cout << "TpcMissAlignment::miniEvent: No simulation mode with" << endl;
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
         if (laser_ray() > 0) {
            return n_sim;
         }
      } else {
         if (beam_point) { // particles produce in one point
            GetPointMuon();
            Int_t passM = passMuoninMF();
            // printf("miniEvent: passMuoninMF=%d\n",passM);
            if (passM > 0) {
               return n_sim;
            }
         } else {
            cout << "TpcMissAlignment::miniEvent: No simulation mode with" << endl;
            cout << "  MFieldMode=" << MFieldMode << "  beam_point=" << beam_point << "  k_par=" << k_par << endl;
            return -2;
         }
      }
   }
   return -1; // number simulation attemps was exeeded
}

//______________________________________________________________________

// laser ray event simulation
Int_t TpcMissAlignment::laser_ray()
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
Int_t TpcMissAlignment::lineClusters(Int_t isec)
{

   Int_t    iRow, iPad, nPads, iPad1, iPad2, nclust;
   Int_t    iRow1, iRow0, iRow2;
   Double_t xloc1, yloc1, xloc2, yloc2, yloc0, ylmin, ylmax, XPads, x1pad, x2pad;
   Double_t UpRowPart, LowRowPart, FirstPadPart, LastPadPart, RowPart, PadPart, spad;
   Double_t hPad, wPad;
   Double_t s1, s2, b1, b2, x11, x12, x21, x22, tng;
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
      if (x11 > x12) {
         Double_t w = x11;
         x11        = x12;
         x12        = w;
      }
      x21 = fR2sec.X() + (b1 - fR2sec.Y()) * tng;
      x22 = fR2sec.X() + (b2 - fR2sec.Y()) * tng;
      if (x21 > x22) {
         Double_t w = x21;
         x21        = x22;
         x22        = w;
      }
      if (x11 > x21) x11 = x21;
      if (x12 < x22) x12 = x22;
      // if(Nreconstructed==107&&isec==3)
      // printf("lineClusters: row=%d xloc(%f,%f)  x(%f,%f)\n",
      // irow,xloc1,xloc2, x11,x12);
      if (xloc1 > x12 || xloc2 < x11) continue;

      if (x11 < xloc1) x11 = xloc1;
      if (x12 > xloc2) x12 = xloc2;
      // if(xloc1>xloc2) printf("lineClusters: error\n");
      nPads = fTpcSecGeo->NofPadsInRow(irow);
      XPads = -nPads * fWpad / 2;
      // printf("lineClusters:     xpad(%f,%f) xloc(%f,%f)  x(%f,%f)\n",
      //  XPads,-XPads,xloc1,xloc2,x11,x12);
      if (x12 < XPads || x11 > -XPads) {
         continue;
      }
      // if(Nreconstructed==11&&isec==2)printf("lineClusters: x11=%f-(%f-%f)*%f=%f\n",fR1sec.X(),b1,fR1sec.Y(),tng,x11);
      // if(Nreconstructed==11&&isec==2)printf("lineClusters: x12=%f-(%f-%f)*%f=%f\n",fR1sec.X(),b2,fR1sec.Y(),tng,x12);
      // if(Nreconstructed==11&&isec==2)printf("lineClusters: x21=%f-(%f-%f)*%f=%f\n",fR2sec.X(),b1,fR2sec.Y(),tng,x21);
      // if(Nreconstructed==11&&isec==2)printf("lineClusters: x22=%f-(%f-%f)*%f=%f\n",fR2sec.X(),b2,fR2sec.Y(),tng,x22);
      if (x11 < XPads) {
         iPad1 = 0;
      } else {
         iPad1 = Int_t((x11 - XPads) / fWpad);
      }
      if (x12 > -XPads) {
         iPad2 = nPads - 1;
      } else { //
               //         if(xloc2>0) iPad2++;
         iPad2 = Int_t((x12 - XPads) / fWpad);
      }
      // if(Nreconstructed==11&&isec==2)printf("lineClusters: irow=%d nPads=%d XPads=%f xloc1,2(%f,%f)
      // iPad1,2(%d->%d)\n",irow,nPads,XPads,xloc1,xloc2,iPad1,iPad2);
      for (Int_t ipad = iPad1; ipad <= iPad2; ipad++) {
         s1          = XPads + ipad * fWpad;
         s2          = s1 + fWpad;
         spad        = padANDband(s1, s2, b1, b2);
         SPads[ipad] = spad * RowPart;
         // if(Nreconstructed==107&&isec==3)
         // printf("lineClusters: %d-%d-%d S(%f,%f)=%f\n",
         // isec,irow,ipad,(s1+s2)/2,(b1+b2)/2,SPads[ipad]);
         // if(Nreconstructed==107&&isec==3)
         // printf("lineClusters: %d-%d-%d (x11,x12,x21,x22)=(%f,%f,%f,%f)\n",isec,irow,ipad,x11,x12,x21,x22);
         // if(Nreconstructed==107&&isec==3)printf("lineClusters: %d-%d-%d   S=%f  rP=%f)\n",
         // isec,irow,ipad,spad,RowPart);
      }
      nclust += ClustersInRow(irow, iPad1, iPad2, iRoc, SPads);
   }
   return nclust;
}

//______________________________________________________________________

// find clusters in the row
Int_t TpcMissAlignment::ClustersInRow(Int_t iRow, Int_t iPad1, Int_t iPad2, Int_t iRoc, Double_t *SPad)
{
   // if(Nreconstructed==10569)
   // printf("  ClustersInRow: s=%d %d(%d,%d)  S(%f ... %f)\n",
   // iSector,iRow,iPad1,iPad2,SPad[iPad1],SPad[iPad2]);
   //
   Double_t Signal, hpad;
   Int_t    nPadsInCluster = 0, nClust = 0;
   hpad = fTpcSecGeo->PadHeight(iRow);
   for (Int_t ipad = iPad1; ipad <= iPad2; ipad++) {
      Signal = RandomPadSignal(iRoc, SPad[ipad]);
      // printf("  ClustersInRow: ipad=%d SPad=%f Signal=%f\n",ipad,SPad[ipad],Signal);
      if (Signal == 0) {
         if (nPadsInCluster > 0) {
            if (newCluster(hpad, nPadsInCluster)) nClust++; // if new cluster created
            // printf("  ClustersInRow:1 ipad=%d nPadsInCluster=%d nClust=%d hits=%d\n",
            // ipad,nPadsInCluster,nClust,MpdTpcMiniMC_Nhits);
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
            zPad[nPadsInCluster] = -simHelixZpos(LCxy, ipad);
         nPadsInCluster++;
         // printf("ClustersInRow:2-1 ipad=%d nPadsInCluster=%d LCxy(%f,%f) z=%f\n",
         // ipad,nPadsInCluster,LCxy.X(),LCxy.Y(),zPad[nPadsInCluster-1]);
         if (nPadsInCluster == maxPcluster) {
            if (newCluster(hpad, nPadsInCluster)) nClust++; // if new cluster created
            // printf("ClustersInRow:2 ipad=%d nPadsInCluster=%d nClust=%d hits=%d\n",
            // ipad,nPadsInCluster,nClust,MpdTpcMiniMC_Nhits);
            nPadsInCluster = 0;
         }
      }
   }
   if (nPadsInCluster > 0) {
      if (newCluster(hpad, nPadsInCluster)) nClust++; // if new cluster created
      // printf("  ClustersInRow:3 nPadsInCluster=%d nClust=%d hits=%d\n",
      // nPadsInCluster,nClust,MpdTpcMiniMC_Nhits);
   }
   return nClust;
}

//______________________________________________________________________

// create new cluster
Bool_t TpcMissAlignment::newCluster(Double_t hpad, Int_t nPads)
{
   Double_t SPAD = 0., xClaster = 0., yClaster = 0., zClaster = 0.;
   Double_t cosAB, Zalong, sigma, sigmac, sigmaz;
   TVector3 eef, A, B, BmA, bma, rlClaster, rgClaster, rbClaster;
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
      <rgClaster> - global coordinates of the cluster
      <rlClaster> - local cluster coordinates in uLCS
      <e_3> - unit vector in GSC along the Z-axis of rLCS (to chamber's outside)
      <eef> - unit vector in GSC of the electric field
      cosAB=(e_3u*eef) for the angle between <A> & <B> (pozitive value)
      <A> - vector along E-field from track point to the sector  plane
      <B> - perpendicular from track point to the sector plane
      cosAB=(e_3u*eef) for the angle between <A> & <B> (negative value)
      <BmA>=<B>-<A>  <bma>={glo2loc}*<AmB>
         xl=Xl+<bma>_x  yl=Yl+<bma>_y  zln=--Zalong*cosAB
         <Xg,Yg,Zg>=R0shift+{loc2glo}*<xln,Yln,zln=-Zalong*cosAB>
   */
   if (e_3.Z() > 0) {
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
   Zalong    = TMath::Abs(zClaster);
   cosAB     = e_3u * eef;
   A         = Zalong * eef;
   B         = (Zalong * cosAB) * e_3u;
   BmA       = B - A;
   bma       = glo2loc * BmA;
   rlClaster = TVector3(xClaster + bma.X(), yClaster + bma.Y(), -Zalong * cosAB);
   rgClaster = uR0shift + uloc2glo * rlClaster;
   // printf("newCluster: Zalong=%f cosAB=%f\n",Zalong,cosAB);
   // printf("newCluster: e_3u(%f,%f,%f)\n",e_3u.X(),e_3u.Y(),e_3u.Z());
   // printf("newCluster: uR0shift(%f,%f,%f)\n",uR0shift.X(),uR0shift.Y(),uR0shift.Z());
   // printf("newCluster: rL(%f,%f,%f)\n",rlClaster.X(),rlClaster.Y(),rlClaster.Z());
   // printf("newCluster: rG(%f,%f,%f)\n",rgClaster.X(),rgClaster.Y(),rgClaster.Z());
   if (rgClaster.Z() < fZmin || rgClaster.Z() > fZmax) return false;
   sigmac = TMath::Sqrt((fWpad * fWpad + hpad * hpad) / nPads);
   sigmaz = TMath::Sqrt((Zalong * sig_Z_pos) * (Zalong * sig_Z_pos) / nPads);
   sigma  = TMath::Sqrt(sigmac * sigmac + sigmaz * sigmaz);
   //    MpdTpcMiniMC_WHitsc[MpdTpcMiniMC_Nhits]=sigma;
   //    MpdTpcMiniMC_WHitsz[MpdTpcMiniMC_Nhits]=sigmaz;
   MpdTpcMiniMC_WHits[MpdTpcMiniMC_Nhits] = sigma;
   MpdTpcMiniMC_iSec[MpdTpcMiniMC_Nhits]  = iSector;
   MpdTpcMiniMC_XHits[MpdTpcMiniMC_Nhits] = rgClaster.X();
   MpdTpcMiniMC_YHits[MpdTpcMiniMC_Nhits] = rgClaster.Y();
   MpdTpcMiniMC_ZHits[MpdTpcMiniMC_Nhits] = rgClaster.Z();
   if (MFieldMode) {
      rbClaster                               = G2B * rgClaster; // GCS->BCS
      MpdTpcMiniMC_XBHits[MpdTpcMiniMC_Nhits] = rbClaster.X();
      MpdTpcMiniMC_YBHits[MpdTpcMiniMC_Nhits] = rbClaster.Y();
      MpdTpcMiniMC_ZBHits[MpdTpcMiniMC_Nhits] = rbClaster.Z();
   }
   MpdTpcMiniMC_Nhits++;
   /*printf("%dth newCluster(s=%d) of %dpads rL(%f,%f,%f) rG(%f,%f,%f)\n",
   MpdTpcMiniMC_Nhits,iSector,nPads,rlClaster.X(),rlClaster.Y(),rlClaster.Z(),
   rgClaster.X(),rgClaster.Y(),rgClaster.Z());*/

   return true;
}

//______________________________________________________________________

// simulate z pozition
Double_t TpcMissAlignment::simLineZpos(TVector2 LCxy)
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
void TpcMissAlignment::crossLinePoints()
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
Bool_t TpcMissAlignment::crossCheck(Int_t iSec, TVector3 Tr1g, TVector3 Tr2g)
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
   R0shift        = rR0_A[isec];
   loc2glo        = r_loc2glo[isec];
   e_3            = TVector3(loc2glo.XZ(), loc2glo.YZ(), loc2glo.ZZ());
   Double_t t_1   = (e_3 * (R0shift - Tr1g)) / (e_3 * eEf); // to 1st sector plane point
   Double_t t_2   = (e_3 * (R0shift - Tr2g)) / (e_3 * eEf); // to 2nd sector plane point
   TVector3 xg_l1 = Tr1g + t_1 * eEf;                       // Point 1 on the sector plane (GCS)
   TVector3 xg_l2 = Tr2g + t_2 * eEf;                       // Point 2 on the sector plane (GCS)
   G0shift        = rG0_A[isec];
   glo2loc        = r_glo2loc[isec];
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
   if (Yes) { // find the band from the direct line track
      fR1_c    = Tr1g;
      fR2_c    = Tr2g;
      uR0shift = uR0_A[isec];
      uloc2glo = u_loc2glo[isec];
      e_3u     = TVector3(u_loc2glo[isec].XZ(), u_loc2glo[isec].YZ(), u_loc2glo[isec].ZZ());
      // printf("crossCheck: uR0shift(%f,%f,%f)\n",uR0shift.X(),uR0shift.Y(),uR0shift.Z());
   }
   // printf("crossCheck: sector=%d  Yes=%d\n",iSec,Yes);
   return Yes;
}

//______________________________________________________________________
// simulate theta angle for cosmic muon
Double_t TpcMissAlignment::theta_cosm()
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
Double_t TpcMissAlignment::theta_beam()
{
   Double_t q = fRandom->Rndm();
   return TMath::ASin(q);
}

//______________________________________________________________________
// simulate pad signal
Double_t TpcMissAlignment::RandomPadSignal(Int_t iRoc, Double_t SPad)
{
   Double_t sigma = SPad * sig_padS[iRoc];
   if (SPad == 0) return 0.;
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
Double_t TpcMissAlignment::padANDband(Double_t s1, Double_t s2, Double_t b1, Double_t b2)
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
Double_t TpcMissAlignment::sTriangle(Int_t line, Double_t s, Double_t b)
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
Double_t TpcMissAlignment::sTrapezoid(Int_t line, Int_t axes, Double_t s, Double_t b1, Double_t b2)
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
Double_t TpcMissAlignment::sParallelogram(Int_t axes, Double_t s1, Double_t s2, Double_t b1, Double_t b2)
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

/*/______________________________________________________________________
// simulate cosmic muon position and direction
  void TpcMissAlignment::GetCosmicMuon() {
    Double_t z1,z2,r,theta,phi;
    z1=-fL_svol;
    z2=fL_svol;
    if(chambers==1) {
      z1=0;
    }
    else {
      if(chambers==2) {
        z2=0;
      }
      else {
        if(chambers!=0) {
          cout<<"MpdTpcMiniMC: GetCosmicMuon: no the mode chambers="<<
                chambers<<endl;
          return;
        }
      }
    }
    r=fr_svol*TMath::Sqrt(fRandom->Rndm());
    phi=TMath::TwoPi()*fRandom->Rndm();
    // random point inside TPC
    r0Beam=TVector3(r*TMath::Cos(phi),r*TMath::Sin(phi),
                    z1+(z2-z1)*fRandom->Rndm());
    if(k_par==2) {
      theta=theta_cosm();
    }
    else {
      theta=theta_beam();
    }
    phi=TMath::TwoPi()*fRandom->Rndm();
    r=TMath::Sin(theta),
    // random muon momentum
    pLas=TVector3(r*TMath::Sin(phi),-TMath::Cos(theta),r*TMath::Cos(phi));
//    pLas=TVector3(-TMath::Cos(theta),r*TMath::Sin(phi),r*TMath::Cos(phi));
    trackPin[0]=r0Beam.X();
    trackPin[1]=r0Beam.Y();
    trackPin[2]=r0Beam.Z();
    trackPin[3]=TMath::ACos(pLas.Z());
    trackPin[4]=TMath::ATan2(pLas.Y(),pLas.X());
    pLas=(p_min+(p_max-p_min)*fRandom->Rndm())*pLas;
 }*/

//______________________________________________________________________
// simulate cosmic muon position and direction
void TpcMissAlignment::GetCosmicMuon()
{
   Double_t x, x0, x02, x1, x2, y, z, z0, z1, z2, zL, zR, t0, a, b, c, d, theta, phi;
   x1 = -3 * fr_svol;
   x2 = 3 * fr_svol;
   z1 = -3 * fL_svol;
   z2 = 3 * fL_svol;
   zL = -fr_svol;
   zR = fL_svol;
   if (chambers == 1) {
      z1 = -2 * fL_svol;
      zL = 0;
   } else {
      if (chambers == 2) {
         z2 = 2 * fL_svol;
         zR = 0;
      } else {
         if (chambers != 0) {
            cout << "MpdTpcMiniMC: GetCosmicMuon: no the mode chambers=" << chambers << endl;
            return;
         }
      }
   }
   // random theta angle
   if (k_par == 2) {
      theta = theta_cosm();
   } else {
      theta = theta_beam();
   }
   // random phi angle
   phi = TMath::TwoPi() * fRandom->Rndm();
   // random muon momentum
   pLas = TVector3(TMath::Sin(theta) * TMath::Cos(phi), TMath::Cos(theta), TMath::Sin(theta) * TMath::Sin(phi));
   while (true) {
      // random point on the plane y=0
      x0 = x1 + (x2 - x1) * fRandom->Rndm();
      z0 = z1 + (z2 - z1) * fRandom->Rndm();
      if (TMath::Abs(pLas.Z()) < 1e-8 // vertical muon
          && z0 >= zL && z0 <= zR && x0 >= -fr_svol && x0 <= fr_svol)
         break; // inside camera
      t0 = (zL - z0) / pLas.Z();
      x  = x0 + pLas.X() * t0;
      y  = pLas.Y() * t0;
      if (TMath::Sqrt(x * x * y * y) <= fr_svol) break; // croses the left side
      t0 = (zR - z0) / pLas.Z();
      x  = x0 + pLas.X() * t0;
      y  = pLas.Y() * t0;
      if (TMath::Sqrt(x * x * y * y) <= fr_svol) break; // croses the right side
      a   = TMath::Sin(theta) * TMath::Sin(theta);      // at^2+2bt+c=0
      x02 = x0 * x0;
      b   = x0 * pLas.X();
      c   = x02 - fr_svol * fr_svol;
      d   = b * b - c;
      if (d < 0) continue; // muon dosn't cross the camera cylinder
      t0 = (-b - TMath::Sqrt(d)) / a;
      z  = z0 + pLas.Z() * t0;
      if (zL <= z && z <= zR) break;
      t0 = (-b + TMath::Sqrt(d)) / a;
      z  = z0 + pLas.Z() * t0;
      if (zL <= z && z <= zR) break;
   }
   r0Beam      = TVector3(x0, 0, z0);
   trackPin[0] = r0Beam.X();
   trackPin[1] = r0Beam.Y();
   trackPin[2] = r0Beam.Z();
   trackPin[3] = TMath::ACos(pLas.Z());
   trackPin[4] = TMath::ATan2(pLas.Y(), pLas.X());
   pLas        = (p_min + (p_max - p_min) * fRandom->Rndm()) * pLas;
}

//______________________________________________________________________
// Get the i-th track hit(cluster)
void TpcMissAlignment::GetHit(Int_t i, double &x, double &y, double &z, double &w)
{
   x = MpdTpcMiniMC_XHits[i];
   y = MpdTpcMiniMC_YHits[i];
   z = MpdTpcMiniMC_ZHits[i];
   w = MpdTpcMiniMC_WHits[i];
}

//______________________________________________________________________
// Chi2 for the clusters fit by the direct line
double TpcMissAlignment_Chi2Line(const double *par)
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
   // printf("Chi2Line: chi2/nh=%f (%f %f %f %f %f)\n",chisq/MpdTpcMiniMC_Nhits,
   // par[0],par[1],par[2],180/TMath::Pi()*par[3],180/TMath::Pi()*par[4]);
   return chisq;
}

//______________________________________________________________________
// find the fit by the direct line
Int_t TpcMissAlignment::fitLine()
{
   // Choose method upon creation between:
   // kMigrad, kSimplex, kCombined,
   // kScan, kFumili
   static ROOT::Minuit2::Minuit2Minimizer minFit(ROOT::Minuit2::kMigrad);

   ROOT::Math::Functor fcn(&TpcMissAlignment_Chi2Line, 5);
   double              step[5] = {0.1, 0.1, 0.1, 0.01, 0.01};
   double              limL[5] = {-133, -133, -170, 0, -pi};
   double              limU[5] = {133, 133, 170, pi, pi};
   double              par[5];
   Int_t               imin, imax, iminx, imaxx, iminy, imaxy, iminz, imaxz;
   Double_t            xmn = 200, xmx = -200, ymn = 200, ymx = -200, zmn = 200, zmx = -200, X = 0, Y = 0, Z = 0;
   Double_t            ptet, pphi;
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

   // if(Nreconstructed==48)  for(Int_t i=0;i<MpdTpcMiniMC_Nhits;i++) {
   // printf("fitLine: e=%d i=%d   r(%f,%f,%f) sector=%d\n",Nreconstructed,i,
   // MpdTpcMiniMC_XHits[i],MpdTpcMiniMC_YHits[i],MpdTpcMiniMC_ZHits[i],MpdTpcMiniMC_iSec[i]);
   // }
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
   ptet        = par[3];
   pphi        = par[4];
   trackPin[5] = 0;
   // printf("par(%f,%f,%f,%f,%f)\n",par[0],par[1],par[2],r2d*par[3],r2d*par[4]);
   for (Int_t i = 0; i < 5; i++) {
      if (i < 3) {
         limL[i] = par[i] - 15;
         limU[i] = par[i] + 15;
      } else {
         limL[i] = par[i] - 0.05;
         limU[i] = par[i] + 0.05;
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
   const double *xs;
   xs = minFit.X();
   if (minFit.MinValue() / MpdTpcMiniMC_Nhits > EvalpHit) {
      //        minFit.PrintResults();
      /*        printf("fitLine:  event=%d MpdTpcMiniMC_Nhits=%d (%f %f %f %f %f)\n",
                     Nreconstructed,MpdTpcMiniMC_Nhits,trackPin[0],trackPin[1],
                     trackPin[2],r2d*trackPin[3],r2d*trackPin[4]);
              printf("Par: 0(%9.4f %9.4f)",X,xs[0]);
              printf(    " 1(%9.4f %9.4f)",Y,xs[1]);
              printf(    " 2(%9.4f %9.4f)",Z,xs[2]);
              printf(    " 3(%9.4f %9.4f)",r2d*ptet,r2d*pphi);
              printf(    " 4(%9.4f %9.4f)\n",r2d*trackPin[4],r2d*trackPout[4]);
              printf("Lim: 0(%9.4f %9.4f)",limL[0],limU[0]);
              printf(    " 1(%9.4f %9.4f)",limL[1],limU[1]);
              printf(    " 2(%9.4f %9.4f)",limL[2],limU[2]);
              printf(    " 3(%9.4f %9.4f)",r2d*limL[3],r2d*limU[3]);
              printf(    " 4(%9.4f %9.4f)\n",r2d*limL[4],r2d*limU[4]);*/
      return 1;
   }
   trackPout[0] = xs[0];
   trackPout[1] = xs[1];
   trackPout[2] = xs[2];
   trackPout[3] = xs[3];
   trackPout[4] = xs[4];
   trackPout[5] = minFit.NCalls();

   /*      trackPout[0]=trackPin[0];
         trackPout[1]=trackPin[1];
         trackPout[2]=trackPin[2];
         trackPout[3]=trackPin[3];
         trackPout[4]=trackPin[4];*/

   // printf("event=%d Par(%9.4f %9.4f %9.4f %9.4f %9.4f)\n",Nreconstructed,
   // trackPout[0],trackPout[1],trackPout[2],r2d*trackPout[3],r2d*trackPout[4]);
   // printf("Par: 0(%9.4f %9.4f) 1(%9.4f %9.4f) 2(%9.4f %9.4f) 3(%9.4f %9.4f) 4(%9.4f %9.4f)\n",
   // trackPin[0],trackPout[0],trackPin[1],trackPout[1],trackPin[2],trackPout[2],r2d*trackPin[3],r2d*trackPout[3],r2d*trackPin[4],r2d*trackPout[4]);
   for (Int_t k = 0; k < 6; k++) trackDiff[k] = trackPout[k] - trackPin[k];
   return 0;
   //    }
   //    else return 1;
}

//______________________________________________________________________
// put chi2 values
void TpcMissAlignment::putChi2L()
{

   //    Double_t chisq2 = 0;
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

      // printf("putChi2: MpdTpcMiniMC_iSec[%d]=%d\n",i,MpdTpcMiniMC_iSec[i]);
      fNhitSec[MpdTpcMiniMC_iSec[i]]++;
      fSecChi2[MpdTpcMiniMC_iSec[i]] += delta * delta;
   }
   Chi2Fit = chisq / MpdTpcMiniMC_Nhits;
}

//  *************************************************************************
//  ************************** magnetic field part **************************
//  *************************************************************************

//______________________________________________________________________
// Simulate muon parameters at a magnetic field
void TpcMissAlignment::GetPointMuon()
{
   /// p_min!=0 => p_min is minimum radius  k_par is muon charge
   ///             charge is random (+-1)
   /// p_min==0 => pLas - muon momentun, k_par is muon charge

   const Double_t R_min = p_min;
   const Double_t R_max = p_max;
   Double_t       pp_min, pp_max, theta1, theta2, theta, phi, rc2B, rcB, phi_s, t, t1, t2, arg;
   if (p_min != 0) {
      q_muon   = fRandom->Rndm() > 0.5 ? 1 : -1;
      q_sign   = q_muon / TMath::Abs(q_muon); // +-1
      c_p      = 333.56 / (TMath::Abs(q_muon) * bMag);
      c_p_m1   = 1 / c_p;
      signRotB = (q_sign > 0) ? -1 : 1; // rotation sign for helix
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
      Double_t pt_min = c_p_m1 * R_min;
      //    Double_t pt_max=c_p_m1*0.8660254038*fr_svol;      // 0.8660254038=sqrt(3)/2
      Double_t pt_max = c_p_m1 * R_max;               // 0.8660254038=sqrt(3)/2
      if (pt_max == 0) pt_max = c_p_m1 * 3 * fr_svol; // 0.8660254038=sqrt(3)/2
      pTMuonB = pt_min + (pt_max - pt_min) * fRandom->Rndm();
      // cout<<"----------------------------------------pTMuonB="<<pTMuonB<<endl;
      rH     = c_p * pTMuonB; // base helix radius
      rHa[0] = rH - dh_track; // inner helix radius
      rHa[1] = rH + dh_track; // outer helix radius
      rHa[2] = rH;            //  helix radius
      pp_min = 0;

      if (rHa[1] < 0.5 * fr_svol) {
         t_end  = twopi;
         pp_max = fL_svol / twopi;
      } else {
         t_end  = 2 * TMath::ASin(0.5 * fr_svol / rHa[1]);
         pp_max = fL_svol / t_end;
      }

      helix_Pz = pp_min + (pp_max - pp_min) * fRandom->Rndm();
      // helix_Pz=(pp_min+pp_max-pp_min)/2;
      phi = TMath::TwoPi() * fRandom->Rndm();
      // phi=0;
      pMuonB = TVector3(pTMuonB * TMath::Cos(phi), pTMuonB * TMath::Sin(phi), helix_Pz);
      pLas   = B2G * pMuonB;
   } else {
      q_sign = q_muon = k_par;
      signRotB        = (q_sign > 0) ? -1 : 1; // rotation sign for helix
      c_p             = 333.56 / (TMath::Abs(q_muon) * bMag);
      c_p_m1          = 1 / c_p;
      pMuonB          = G2B * pLas;
      pTMuonB         = TMath::Sqrt(pMuonB.X() * pMuonB.X() + pMuonB.Y() * pMuonB.Y());
      rH              = c_p * pTMuonB; // base helix radius
      rHa[0]          = rH - dh_track; // inner helix radius
      rHa[1]          = rH + dh_track; // outer helix radius
      rHa[2]          = rH;            //  helix radius
      if (rHa[1] < 0.5 * fr_svol)
         t_end = twopi;
      else
         t_end = 2 * TMath::ASin(0.5 * fr_svol / rHa[1]);
      helix_Pz = rH * pMuonB.Z() / pTMuonB;
   }

   // Helix parameters
   r0BeamB = G2B * r0Beam; // production point in BCS
   xcB     = r0BeamB.X() + rH * pMuonB.Y() / (q_muon * pTMuonB);
   ycB     = r0BeamB.Y() - rH * pMuonB.X() / (q_muon * pTMuonB);
   rc2B    = xcB * xcB + ycB * ycB;
   rcB     = TMath::Sqrt(rc2B);
   rc0_B   = TVector3(xcB, ycB, 0); // point on the cilinder axis in BSC
   rc0_G   = B2G * rc0_B;           // point on the cilinder axis in GSC
   phi_0B  = TMath::ATan2(-ycB, -xcB);
   // printf("pMuon(%f,%f,%f) rH=%f\n",pLas.X(),pLas.Y(),pLas.Z(),rH*q_sign);
   // printf("pMuonB(%f,%f,%f) rH=%f\n",pMuonB.X(),pMuonB.Y(),pMuonB.Z(),rH);
   // printf("(parHelix(%6.1f,%6.1f,%6.1f,%6.1f,%8.1f) %f\n",xcB,ycB,q_sign*rH,
   // phi_0B*r2d,helix_Pz,t_end*r2d);
   // if(simEvt==5) exit(0);
   trackPin[0] = xcB;
   trackPin[1] = ycB;
   trackPin[2] = r0BeamB.Z();
   trackPin[3] = q_muon * rH;
   trackPin[4] = phi_0B;
   trackPin[5] = helix_Pz;
}

//______________________________________________________________________
// simulate pads signals for a muon in the magnetic field
Int_t TpcMissAlignment::passMuoninMF()
{
   // printf("TpcMissAlignment::passMuoninMF()  t_end=%f================\n",t_end);
   Bool_t   InPad;
   Int_t    nclust, irow, irow0, ipad, ipad0, ip1, ip2, nPads, iSector0, sectorI, sec1;
   Double_t t, t0, ts, t0_end, dt, step = 0.1, dtd2, tc, xt, xt0, xc, x1, x2, yt, yt0, yc, y1, y2, hPad, wPad;
   MpdTpcMiniMC_Nhits = 0;
   dt                 = step / rH; // step along phi
   dtd2               = 0.5 * dt;
   //  t0_end=t_end;
   sectorI = iSector = sector0 = coordinates(0, rHa[2]);
   nsPass                      = 0;
   //  InPad=false;              // particle did not cross any sensitive sector area
   //  while(t_S[0]<t_end||t_S[1]<t_end||t_S[2]<t_end) {  // loop along full helix
   // printf("while(t<=t_end)   t=%f  t0=%f t_end=%f\n",t*r2d,t_end*r2d,t_end*r2d);
   for (Int_t s = 0; s < 24; s++) {
      secPass[s] = false;
      for (Int_t i = 0; i < mxp; i++)
         for (Int_t j = 0; j < mxr; j++)
            for (Int_t k = 0; k < 3; k++) {
               SindexPad[s][i][j][k] = 0;
               Sxc1[s][i][j][k]      = 0;
               Sxc2[s][i][j][k]      = 0;
               Syc1[s][i][j][k]      = 0;
               Syc2[s][i][j][k]      = 0;
               StPad1[s][i][j][k]    = 0;
               StPad2[s][i][j][k]    = 0;
            }
   }
   for (Int_t s = 0; s < 12; s++) {
      // sector aligment
      iSector      = s;
      R0shift      = rR0_A[s];
      loc2glo      = r_loc2glo[s];
      e_3          = TVector3(loc2glo.XZ(), loc2glo.YZ(), loc2glo.ZZ());
      G0shift      = rG0_A[s];
      glo2loc      = r_glo2loc[s];
      uR0shift     = uR0_A[s];
      uloc2glo     = u_loc2glo[s];
      e_3u         = TVector3(u_loc2glo[s].XZ(), u_loc2glo[s].YZ(), u_loc2glo[s].ZZ());
      TVector3 hcB = TVector3(xcB, ycB, 0);
      TVector3 hcL = glo2loc * hcB;
      xcL          = hcL.X();
      ycL          = hcL.Y();

      for (Int_t l = 2; l >= 0; l--) { // 3 trips along one sector
                                       // l=2,1,0 own particle, inner and outer circles
                                       // printf("for l=%d  %f<%f\n",l,t_S[l],t_end*r2d);
         //     if(t_S[l]<t_end) {    // circle trip is not finished
         t    = dt;
         irow = -1;
         ipad = -1;
         while (t <= t_end) {
            t0   = t;
            rLt0 = rLt;
            // if(l==2)printf("l=%d t=%f  r0(%f,%f) sec0=%d\n",l,t*r2d,rLt0.X(),rLt0.Y(),s);
            coordinates(s, t, rHa[l]);
            if (rLt.Y() > -0.5) break;
            t += dt;
         }
         while (t <= t_end) {
            t0   = t;
            rLt0 = rLt;
            t += dt;
            coordinates(s, t, rHa[l]);
            // if(simEvt==4)
            // printf("   1 l=%d t=%f  r0(%f,%f) sec[%d,%d]\n",l,t*r2d,rLt.X(),rLt.Y(),sector0,iSector);
            //           if(iSector==sector0) {                         // step inside a sector

            /*if(iSector==11)
            printf("   2 l=%d t=%f  r(%f,%f) sec=%d r0p0i[%d,%d] rpi[%d,%d]\n",
            l,t*r2d,rLt.X(),rLt.Y(),iSector,irow0,ipad0,irow,ipad);*/
            xt0   = rLt0.X();
            yt0   = rLt0.Y();
            xt    = rLt.X();
            yt    = rLt.Y();
            ipad0 = ipad;
            irow0 = irow;
            if (fTpcSecGeo->IsItPad(xt, yt, irow, ipad, x1, x2, y1, y2)) {
               secPass[s] = true;
               // the point is inside the pad(irow,ipad)
               //        sectorI=iSector;
               x1 = xminPad[ipad][irow];
               y1 = yminPad[ipad][irow];
               x2 = xmaxPad[ipad][irow];
               y2 = ymaxPad[ipad][irow];

               if (ipad != ipad0 || irow != irow0) {     // just a NEXT pad
                                                         /*if(3==3)
                                                         printf(" NewPad l=%d t=%f sec=%d r0p0[%d,%d] rp[%d,%d] rLt0(%f,%f) rLt(%f,%f)\n",
                                                         l,t*r2d,s,irow0,ipad0,irow,ipad,xt0,yt0,xt,yt);*/
                                                         // check the bottom border crossing
                                                         // ********************************
                  if (rLt0.Y() <= y1 && y1 <= rLt.Y()) { // track comes from below row
                     tc = (y1 - yt0) / (yt - yt0);
                     xc = xt0 + (xt - xt0) * tc;
                     if (x1 <= xc && xc <= x2) {
                        // track crosses the lower pad bound
                        Sxc1[s][ipad][irow][l]      = xc;
                        Syc1[s][ipad][irow][l]      = y1;
                        StPad1[s][ipad][irow][l]    = t0 + dt * tc;
                        SindexPad[s][ipad][irow][l] = 1; // bottom border cross index

                        if (irow0 >= 0) { // the previous point was inside other pad
                           if (xminPad[ipad0][irow0] <= xc && xc <= xmaxPad[ipad0][irow0]) {
                              Sxc2[s][ipad0][irow0][l]   = xc;
                              Syc2[s][ipad0][irow0][l]   = y1;
                              StPad2[s][ipad0][irow0][l] = t0 + dt * tc;
                              SindexPad[s][ipad0][irow0][l] += 100; // top exit
                           } else {
                              if (xc <= xminPad[ipad0][irow0]) { // came from r-b pad
                                 tc                         = (x2 - xt0) / (xt - xt0);
                                 yc                         = yt0 + (yt - yt0) * tc;
                                 Sxc2[s][ipad0][irow0][l]   = x1;
                                 Syc2[s][ipad0][irow0][l]   = yc;
                                 StPad2[s][ipad0][irow0][l] = t0 + dt * tc;
                                 SindexPad[s][ipad0][irow0][l] += 10; // l-exit
                              } else {                                // came from r-b pad
                                 tc                         = (x1 - xt0) / (xt - xt0);
                                 yc                         = yt0 + (yt - yt0) * tc;
                                 Sxc2[s][ipad0][irow0][l]   = x2;
                                 Syc2[s][ipad0][irow0][l]   = yc;
                                 StPad2[s][ipad0][irow0][l] = t0 + dt * tc;
                                 SindexPad[s][ipad0][irow0][l] += 1000; // r-exit
                              }
                           }
                        } // end of the bottom entry from the previous pad
                     } else {
                        if (xc <= x1) {
                           tc                          = (x1 - xt0) / (xt - xt0);
                           yc                          = yt0 + (yt - yt0) * tc;
                           Sxc1[s][ipad][irow][l]      = x1;
                           Syc1[s][ipad][irow][l]      = yc;
                           StPad1[s][ipad][irow][l]    = t - dtd2;
                           SindexPad[s][ipad][irow][l] = 10; // left border cross index
                           if (irow0 >= 0) {                 // the previous point was inside other pad
                              Sxc2[s][ipad0][irow0][l]   = xc;
                              Syc2[s][ipad0][irow0][l]   = y1;
                              StPad2[s][ipad0][irow0][l] = t0 + dt * tc;
                              SindexPad[s][ipad0][irow0][l] += 100; // t-exit
                           }
                        } else {
                           tc                          = (x2 - xt0) / (xt - xt0);
                           yc                          = yt0 + (yt - yt0) * tc;
                           Sxc1[s][ipad][irow][l]      = x2;
                           Syc1[s][ipad][irow][l]      = yc;
                           StPad2[s][ipad][irow][l]    = t - dtd2;
                           SindexPad[s][ipad][irow][l] = 1000; // right border cross index
                           if (irow0 >= 0) {                   // the previous point was inside other pad
                              Sxc2[s][ipad0][irow0][l]   = xc;
                              Syc2[s][ipad0][irow0][l]   = y1;
                              StPad2[s][ipad0][irow0][l] = t0 + dt * tc;
                              SindexPad[s][ipad0][irow0][l] += 100; // t-exit
                           }
                        }
                     }
                  }
                  // ********************************
                  else {
                     if (rLt.Y() <= y2 && y2 <= rLt0.Y()) {
                        // track comes from above row
                        tc = (y2 - yt0) / (yt - yt0);
                        xc = xt0 + (xt - xt0) * tc;
                        if (x1 <= xc && xc <= x2) {
                           // track crosses the top pad bound
                           Sxc1[s][ipad][irow][l]      = xc;
                           Syc1[s][ipad][irow][l]      = y2;
                           StPad1[s][ipad][irow][l]    = t0 + dt * tc;
                           SindexPad[s][ipad][irow][l] = 100; // top border cross index
                           if (irow0 >= 0) {                  // the previous point was inside other pad
                              if (xminPad[ipad0][irow0] <= xc && xc <= xmaxPad[ipad0][irow0]) {
                                 Sxc2[s][ipad0][irow0][l]   = xc;
                                 Syc2[s][ipad0][irow0][l]   = y2;
                                 StPad2[s][ipad0][irow0][l] = t0 + dt * tc;
                                 SindexPad[s][ipad0][irow0][l] += 1; // bottom exit
                              } else {
                                 if (xc <= xminPad[ipad0][irow0]) { // came from l-t pad
                                    tc                         = (x2 - xt0) / (xt - xt0);
                                    yc                         = yt0 + (yt - yt0) * tc;
                                    Sxc2[s][ipad0][irow0][l]   = x1;
                                    Syc2[s][ipad0][irow0][l]   = yc;
                                    StPad2[s][ipad0][irow0][l] = t0 + dt * tc;
                                    SindexPad[s][ipad0][irow0][l] += 10; // l-exit
                                 } else {                                // came from r-b pad
                                    tc                         = (x1 - xt0) / (xt - xt0);
                                    yc                         = yt0 + (yt - yt0) * tc;
                                    Sxc2[s][ipad0][irow0][l]   = x2;
                                    Syc2[s][ipad0][irow0][l]   = yc;
                                    StPad2[s][ipad0][irow0][l] = t0 + dt * tc;
                                    SindexPad[s][ipad0][irow0][l] += 1000; // r-exit
                                 }
                              }
                           } // end of the top entry from the previous pad
                        } else {
                           if (xc <= x1) {
                              tc                          = (x1 - xt0) / (xt - xt0);
                              yc                          = yt0 + (yt - yt0) * tc;
                              Sxc1[s][ipad][irow][l]      = x1;
                              Syc1[s][ipad][irow][l]      = yc;
                              StPad1[s][ipad][irow][l]    = t0 + dt * tc;
                              SindexPad[s][ipad][irow][l] = 10; // left border cross index
                              if (irow0 >= 0) {                 // the previous point was inside other pad
                                 Sxc2[s][ipad0][irow0][l]   = xc;
                                 Syc2[s][ipad0][irow0][l]   = y1;
                                 StPad2[s][ipad0][irow0][l] = t0 + dt * tc;
                                 SindexPad[s][ipad0][irow0][l] += 1; // bottom exit
                              }
                           } else {
                              tc                          = (x2 - xt0) / (xt - xt0);
                              yc                          = yt0 + (yt - yt0) * tc;
                              Sxc1[s][ipad][irow][l]      = x2;
                              Syc1[s][ipad][irow][l]      = yc;
                              StPad1[s][ipad][irow][l]    = t0 + dt * tc;
                              SindexPad[s][ipad][irow][l] = 1000; // right border cross index
                              if (irow0 >= 0) {                   // the previous point was inside other pad
                                 Sxc2[s][ipad0][irow0][l]   = xc;
                                 Syc2[s][ipad0][irow0][l]   = y1;
                                 StPad2[s][ipad0][irow0][l] = t0 + dt * tc;
                                 SindexPad[s][ipad0][irow0][l] += 1; // bottom exit
                              }
                           }
                        }
                     } else {
                        // ********************************
                        if (rLt0.X() <= x1) {
                           tc                          = (x1 - xt0) / (xt - xt0);
                           yc                          = yt0 + (yt - yt0) * tc;
                           Sxc1[s][ipad][irow][l]      = x1;
                           Syc1[s][ipad][irow][l]      = yc;
                           StPad1[s][ipad][irow][l]    = t0 + dt * tc;
                           SindexPad[s][ipad][irow][l] = 10; // left border cross index
                           if (irow0 >= 0) {                 // the previous point was inside other pad
                              Sxc2[s][ipad0][irow0][l]   = x1;
                              Syc2[s][ipad0][irow0][l]   = yc;
                              StPad2[s][ipad0][irow0][l] = t0 + dt * tc;
                              SindexPad[s][ipad0][irow0][l] += 1000; // right exit
                           }
                        } else {
                           tc                          = (x2 - xt0) / (xt - xt0);
                           yc                          = yt0 + (yt - yt0) * tc;
                           Sxc1[s][ipad][irow][l]      = x2;
                           Syc1[s][ipad][irow][l]      = yc;
                           StPad1[s][ipad][irow][l]    = t0 + dt * tc;
                           SindexPad[s][ipad][irow][l] = 1000; // right border cross index
                           if (irow0 >= 0) {                   // the previous point was inside other pad
                              Sxc2[s][ipad0][irow0][l]   = x2;
                              Syc2[s][ipad0][irow0][l]   = yc;
                              StPad2[s][ipad0][irow0][l] = t0 + dt * tc;
                              SindexPad[s][ipad0][irow0][l] += 10; // left exit
                           }
                        }
                     }
                  } // ********************************
               }    // check of ipad0!=ipad
                    // if(iSector==10)
               // printf("   3 l=%d t=%f  r(%f,%f) sec=%d r0p0[%d,%d] rpi[%d,%d  %d]\n",
               // l,t*r2d,rLt.X(),rLt.Y(),iSector,irow0,ipad0,irow,ipad,indexPad[ipad][irow][l]);
            }      // check of IsItPad
            else { // track is out of pad
               /*if(t*r2d>32.1)printf("-> l=%d sec=%d t=%f rt0(%f,%f) rt(%f,%f) r,p[%d,%d-%d,%d]\n",
                   l,iSector,t*r2d,xt0,yt0,xt,yt,irow0,ipad0,irow,ipad);*/
               if (ipad >= 0) { // exit from the pad(ipad0,irow0)
                  // find the border of the exit
                  // x1,y1,x2,y2 for pad(ipad0,irow0)
                  x1 = xminPad[ipad0][irow0];
                  y1 = yminPad[ipad0][irow0];
                  x2 = xmaxPad[ipad0][irow0];
                  y2 = ymaxPad[ipad0][irow0];
                  // check the bottom border crossing
                  if (rLt.Y() <= y1 && y1 <= rLt0.Y()) {
                     tc = (y1 - yt0) / (yt - yt0);
                     xc = xt0 + (xt - xt0) * tc;
                     if (x1 <= xc && xc <= x2) {
                        // exit from bottom
                        Sxc2[s][ipad0][irow0][l]   = xc;
                        Syc2[s][ipad0][irow0][l]   = y1;
                        StPad2[s][ipad0][irow0][l] = t0 + dt * tc;
                        SindexPad[s][ipad0][irow0][l] += 1; // bottom border cross index
                     } else {
                        if (xc <= xminPad[ipad0][irow0]) { // came from l-b pad
                           tc                         = (x1 - xt0) / (xt - xt0);
                           yc                         = yt0 + (yt - yt0) * tc;
                           Sxc2[s][ipad0][irow0][l]   = x1;
                           Syc2[s][ipad0][irow0][l]   = yc;
                           StPad2[s][ipad0][irow0][l] = t0 + dt * tc;
                           SindexPad[s][ipad0][irow0][l] += 10; // l-exit
                        } else {                                // came from r-b pad
                           tc                         = (x2 - xt0) / (xt - xt0);
                           yc                         = yt0 + (yt - yt0) * tc;
                           Sxc2[s][ipad0][irow0][l]   = x2;
                           Syc2[s][ipad0][irow0][l]   = yc;
                           StPad2[s][ipad0][irow0][l] = t0 + dt * tc;
                           SindexPad[s][ipad0][irow0][l] += 1000; // r-exit
                        }
                     }
                  } // end of the bottom exit from the pad
                  else {
                     if (rLt0.Y() <= y2 && y2 <= rLt.Y()) {
                        tc = (y2 - yt0) / (yt - yt0);
                        xc = xt0 + (xt - xt0) * tc;
                        if (x1 <= xc && xc <= x2) { // exit from the top
                           Sxc2[s][ipad0][irow0][l]   = xc;
                           Syc2[s][ipad0][irow0][l]   = y2;
                           StPad2[s][ipad0][irow0][l] = t0 + dt * tc;
                           SindexPad[s][ipad0][irow0][l] += 100;
                        } else {
                           if (xc <= xminPad[ipad0][irow0]) { // exit on the left
                              tc                         = (x1 - xt0) / (xt - xt0);
                              yc                         = yt0 + (yt - yt0) * tc;
                              Sxc2[s][ipad0][irow0][l]   = x1;
                              Syc2[s][ipad0][irow0][l]   = yc;
                              StPad2[s][ipad0][irow0][l] = t0 + dt * tc;
                              SindexPad[s][ipad0][irow0][l] += 10;
                           } else { // exit on the right
                              tc                         = (x2 - xt0) / (xt - xt0);
                              yc                         = yt0 + (yt - yt0) * tc;
                              Sxc2[s][ipad0][irow0][l]   = x2;
                              Syc2[s][ipad0][irow0][l]   = yc;
                              StPad2[s][ipad0][irow0][l] = t0 + dt * tc;
                              SindexPad[s][ipad0][irow0][l] += 1000;
                           }
                        }
                     } // end of the top exit from the pad
                     else {
                        if (rLt.X() <= x1 && x1 <= rLt0.X()) {
                           tc                         = (x1 - xt0) / (xt - xt0);
                           yc                         = yt0 + (yt - yt0) * tc;
                           Sxc2[s][ipad0][irow0][l]   = x1;
                           Syc2[s][ipad0][irow0][l]   = yc;
                           StPad2[s][ipad0][irow0][l] = t0 + dt * tc;
                           SindexPad[s][ipad0][irow0][l] += 10; // left border cross index
                        } else {
                           tc                         = (x2 - xt0) / (xt - xt0);
                           yc                         = yt0 + (yt - yt0) * tc;
                           Sxc2[s][ipad0][irow0][l]   = x2;
                           Syc2[s][ipad0][irow0][l]   = yc;
                           StPad2[s][ipad0][irow0][l] = t0 + dt * tc;
                           SindexPad[s][ipad0][irow0][l] += 1000; // right border index
                        }
                     }
                  }
                  ipad = -1;
                  irow = -1;
               }
            } // end of "not a pad"
         }    // end of while(t<t_end)
      }       // end of l-loop
              //    for(Int_t s=0;s<24;s++) {               // check sectors for a hit
      if (secPass[s]) {
         for (Int_t i = 0; i < mxp; i++)
            for (Int_t j = 0; j < mxr; j++)
               for (Int_t k = 0; k < 3; k++) {
                  indexPad[i][j][k] = SindexPad[s][i][j][k];
                  /*if(k<2)*/ {
                     xc1[i][j][k]   = Sxc1[s][i][j][k];
                     xc2[i][j][k]   = Sxc2[s][i][j][k];
                     yc1[i][j][k]   = Syc1[s][i][j][k];
                     yc2[i][j][k]   = Syc2[s][i][j][k];
                     tPad1[i][j][k] = StPad1[s][i][j][k];
                     tPad2[i][j][k] = StPad2[s][i][j][k];
                  }
               }
         CheckPadsinMF();
         // ---------------------
         if (2 == 1) {
            Int_t nRows    = fTpcSecGeo->NofRows();
            Int_t maxnPads = fTpcSecGeo->NofPadsInRow(nRows - 1);
            printf("         ============================================");
            printf(" sector=%d ============================================\n", s);
            for (Int_t j = nRows - 1; j >= 0; j--) {
               Int_t rPads  = fTpcSecGeo->NofPadsInRow(j);
               Int_t wsPads = (maxnPads - rPads) / 2;
               for (Int_t i = 0; i < wsPads; i++) printf(" ");
               for (Int_t i = 0; i < rPads; i++)
                  if (indexPad[i][j][0] != 0 || indexPad[i][j][1] != 0 || indexPad[i][j][2] != 0)
                     printf("*");
                  else
                     printf("0");
               cout << endl;
               sector1 = s;
            }
         }
         // ---------------------
      }
   } // end of while(t<t_end)
   return 1;
   // printf("passMounInMF: MpdTpcMiniMC_Nhits=%d\n",MpdTpcMiniMC_Nhits);
   // exit(0);
}

//   ***************************************************
//   ***************************************************
//   ***************************************************
//______________________________________________________________________
// find all clusters in a sector from a particle in the magnetic field
void TpcMissAlignment::CheckPadsinMF()
{
   Int_t    nclust, irow, ipad, ipl0[2][2], ipl1[2][2], ip0, ip1, ip2, nPads, nRows, n0, n1, ipw1, ipw2, ipt;
   Double_t hPad, wPad, SPad, sw1, sw2, x1, y1, x2, y2;
   //   ***************************************************
   //   ********       check pads in each row      ********
   //   ***************************************************

   nclust = 0;
   wPad   = fWpad;
   nRows  = fTpcSecGeo->NofRows();
   // printf("==================== Sector=%d ======================\n",iSector);
   for (irow = 0; irow < nRows; irow++) { //  find clusters in the row
      hPad  = fTpcSecGeo->PadHeight(irow);
      SPad  = wPad * hPad;
      nPads = fTpcSecGeo->NofPadsInRow(irow);
      // printf("passMuoninMF: *************** irow=%d ***** nPads=%d************\n",
      // irow,nPads);

      for (Int_t i = 0; i < 2; i++)
         for (Int_t j = 0; j < 2; j++) {
            ipl0[i][j] = -1;
            ipl1[i][j] = -1;
         }
      ip0 = -10;
      n0  = -1;
      for (ipad = 0; ipad < nPads; ipad++) {
         if (indexPad[ipad][irow][0] != 0) {
            // printf("p0=%d index=%d\n",ipad,indexPad[ipad][irow][0]);
            if (ip0 < 0) {
               n0++;
               ip0 = ipl0[n0][0] = ipl0[n0][1] = ipad;
            } else {
               ipl0[n0][1] = ipad;
            }
         } else {
            if (ip0 > 0 && n0 == 0) ip0 = -10;
         }
      }
      ip0 = -10;
      n1  = -1;
      n0++;
      for (ipad = 0; ipad < nPads; ipad++) {
         if (indexPad[ipad][irow][1] != 0) {
            // printf("p1=%d index=%d\n",ipad,indexPad[ipad][irow][1]);
            if (ip0 < 0) {
               n1++;
               ip0 = ipl1[n1][0] = ipl1[n1][1] = ipad;
            } else {
               ipl1[n1][1] = ipad;
            }
         } else {
            if (ip0 > 0 && n1 == 0) ip0 = -10;
         }
      }
      n1++;
      /*printf("n0=%d n1=%d ipl0[%d,%d %d,%d] ipl1[%d,%d %d,%d]]\n",
      n0,n1,ipl0[0][0],ipl0[0][1],ipl0[1][0],ipl0[1][1],
      ipl1[0][0],ipl1[0][1],ipl1[1][0],ipl1[1][1]);*/

      if (n0 == 0 && n1 == 0) {
         ip1 = 2;
         ip2 = 1;
      } else {
         if ((n0 == 1 && n1 == 0) || (n0 == 0 && n1 == 1)) {
            n0orn1is1(n0, n1, irow, nPads, ipl0, ipl1, ip1, ip2);
         } else {
            if (n0 == 1 && n1 == 1) {
               n0is1n1is1(irow, ipl0, ipl1, ip1, ip2);
            } else {
               if (n0 == 2 && n1 == 1) {
                  if (ipl1[0][0] < ipl0[1][0]) { // case: "/../" & "\"
                     n0is1n1is1(irow, ipl0, ipl1, ip1, ipw2);
                     ipl0[0][0] = ipl0[1][0];
                     ipl0[0][1] = ipl0[1][1];
                     n0orn1is1(1, n0, irow, nPads, ipl0, ipl1, ipw1, ip2);
                  } else {
                     n0is1n1is1(irow, ipl0, ipl1, ipw1, ip2);
                     n0orn1is1(1, n0, irow, nPads, ipl0, ipl1, ip1, ipw2);
                  }
               }
            }
         }
      }
      // printf("passMuoninMF: rp[%d,%d-%d]\n",irow,ip1,ip2);
      for (ipad = ip1; ipad <= ip2; ipad++) {
         x1 = xminPad[ipad][irow];
         y1 = yminPad[ipad][irow];
         x2 = xmaxPad[ipad][irow];
         y2 = ymaxPad[ipad][irow];

         /*if(indexPad[ipad][irow][0]!=0||indexPad[ipad][irow][1]!=0
         ||indexPad[ipad][irow][2]!=0)
         printf("CheckPadsinMF: sec=%2d indexPad[%2d,%2d](0,1,2)=(%4d %4d %4d)\n",iSector,
         irow,ipad,indexPad[ipad][irow][0],indexPad[ipad][irow][1],indexPad[ipad][irow][2]);*/

         //   ***************************************************
         //             valid indexes combinations
         // 0000  1001 I 0011  0101 I 0101  0110 I 1010  0110
         // 0000  0101 I 0011  0011 I 0101  1100 I 1010  1010
         // 0000  0011 I 0011  1010 I 1001  0101 I 1010  1100
         // 0000  0110 I 0011  1100 I 1001  1001 I 0110  0110
         // 0000  1010 I 0101  1001 I 1001  0110 I 0110  1010
         // 0000  1100 I 0101  0101 I 1001  1010 I 1100  1100
         //
         //   ***************************************************
         if (indexPad[ipad][irow][0] == 2000) indexPad[ipad][irow][0] = 0;
         if (indexPad[ipad][irow][1] == 2000) indexPad[ipad][irow][1] = 0;
         if (indexPad[ipad][irow][0] == 20) indexPad[ipad][irow][0] = 0;
         if (indexPad[ipad][irow][1] == 20) indexPad[ipad][irow][1] = 0;

         if (indexPad[ipad][irow][0] == 0 && indexPad[ipad][irow][1] == 0) {
            if (indexPad[ipad][irow][2] != 0) { // case a) of padANDband
               SPads[ipad] = wPad * hPad;
            } else {
               SPads[ipad] = 0;
            }
            continue;
         }

         // case of one curve cross (index1==0, index2!=00  )
         //
         // case d) 0-1001
         if (indexPad[ipad][irow][0] == 0 && indexPad[ipad][irow][1] == 1001) {
            if (xc1[ipad][irow][1] < x2)
               SPads[ipad] = 0.5 * (x2 - xc1[ipad][irow][1]) * (yc2[ipad][irow][1] - y1);
            else
               SPads[ipad] = 0.5 * (x2 - xc2[ipad][irow][1]) * (yc1[ipad][irow][1] - y1);
            continue;
         }
         // case b) 1001-0
         if (indexPad[ipad][irow][1] == 0 && indexPad[ipad][irow][0] == 1001) {
            if (xc1[ipad][irow][0] < x2)
               sw1 = 0.5 * (x2 - xc1[ipad][irow][0]) * (yc2[ipad][irow][0] - y1);
            else
               sw1 = 0.5 * (x2 - xc2[ipad][irow][0]) * (yc1[ipad][irow][0] - y1);
            SPads[ipad] = SPad - sw1;
            continue;
         }
         // case c) 0-0101
         if (indexPad[ipad][irow][0] == 0 && indexPad[ipad][irow][1] == 101) {
            Bool_t centerONleft;
            if (yc1[ipad][irow][1] < yc2[ipad][irow][1]) {
               if (xc1[ipad][irow][1] < xc2[ipad][irow][1])
                  centerONleft = false;
               else
                  centerONleft = true;
            } else {
               if (xc1[ipad][irow][1] < xc2[ipad][irow][1])
                  centerONleft = true;
               else
                  centerONleft = false;
            }
            if (centerONleft)
               SPads[ipad] = 0.5 * (xc1[ipad][irow][1] + xc2[ipad][irow][1] - 2 * x1) * hPad;
            else
               SPads[ipad] = 0.5 * (2 * x2 - xc1[ipad][irow][1] - xc2[ipad][irow][1]) * hPad;
            continue;
         }
         // case c) 0101-0
         if (indexPad[ipad][irow][0] == 101 && indexPad[ipad][irow][1] == 0) {
            Bool_t centerONleft;
            if (yc1[ipad][irow][0] < yc2[ipad][irow][0]) {
               if (xc1[ipad][irow][0] < xc2[ipad][irow][0])
                  centerONleft = false;
               else
                  centerONleft = true;
            } else {
               if (xc1[ipad][irow][0] < xc2[ipad][irow][0])
                  centerONleft = true;
               else
                  centerONleft = false;
            }
            if (centerONleft)
               SPads[ipad] = 0.5 * (2 * x2 - xc1[ipad][irow][0] - xc2[ipad][irow][0]) * hPad;
            else
               SPads[ipad] = 0.5 * (xc1[ipad][irow][0] + xc2[ipad][irow][0] - 2 * x1) * hPad;
            continue;
         }
         // case c) 0-0101
         if (indexPad[ipad][irow][0] == 0 && indexPad[ipad][irow][1] == 101) {
            Bool_t centerONleft;
            if (yc1[ipad][irow][1] < yc2[ipad][irow][1]) {
               if (xc1[ipad][irow][1] < xc2[ipad][irow][1])
                  centerONleft = false;
               else
                  centerONleft = true;
            } else {
               if (xc1[ipad][irow][1] < xc2[ipad][irow][1])
                  centerONleft = true;
               else
                  centerONleft = false;
            }
            if (centerONleft)
               SPads[ipad] = 0.5 * (2 * x2 - xc1[ipad][irow][1] - xc2[ipad][irow][1]) * hPad;
            else
               SPads[ipad] = 0.5 * (xc1[ipad][irow][1] + xc2[ipad][irow][1] - 2 * x1) * hPad;
            continue;
         }
         // case b) 0-0011
         if (indexPad[ipad][irow][0] == 0 && indexPad[ipad][irow][1] == 11) {
            if (xc1[ipad][irow][1] > x1)
               sw1 = 0.5 * (xc1[ipad][irow][1] - x1) * (yc2[ipad][irow][1] - y1);
            else
               sw1 = 0.5 * (xc2[ipad][irow][1] - x1) * (yc1[ipad][irow][1] - y1);
            SPads[ipad] = SPad - sw1;
            continue;
         }
         // case d) 0011-0
         if (indexPad[ipad][irow][1] == 0 && indexPad[ipad][irow][0] == 11) {
            if (xc1[ipad][irow][0] > x1)
               SPads[ipad] = 0.5 * (xc1[ipad][irow][0] - x1) * (yc2[ipad][irow][0] - y1);
            else
               SPads[ipad] = 0.5 * (xc2[ipad][irow][0] - x1) * (yc1[ipad][irow][0] - y1);
            continue;
         }
         // case b) 0-0110
         if (indexPad[ipad][irow][0] == 0 && indexPad[ipad][irow][1] == 110) {
            if (xc1[ipad][irow][1] > x1)
               sw1 = 0.5 * (xc1[ipad][irow][1] - x1) * (y2 - yc2[ipad][irow][1]);
            else
               sw1 = 0.5 * (xc2[ipad][irow][1] - x1) * (y2 - yc1[ipad][irow][1]);
            SPads[ipad] = SPad - sw1;
            continue;
         }
         // case d) 0110-0
         if (indexPad[ipad][irow][1] == 0 && indexPad[ipad][irow][0] == 110) {
            if (xc1[ipad][irow][0] > x1)
               SPads[ipad] = 0.5 * (xc1[ipad][irow][0] - x1) * (y2 - yc2[ipad][irow][0]);
            else
               SPads[ipad] = 0.5 * (xc2[ipad][irow][0] - x1) * (y2 - yc1[ipad][irow][0]);
            continue;
         }
         // case c) 0-1010
         if (indexPad[ipad][irow][0] == 0 && indexPad[ipad][irow][1] == 1010) {
            SPads[ipad] = 0.5 * (yc1[ipad][irow][1] + yc2[ipad][irow][1] - 2 * y1) * (x2 - x1);
            continue;
         }
         // case c) 0110-0
         if (indexPad[ipad][irow][1] == 0 && indexPad[ipad][irow][0] == 1010) {
            SPads[ipad] = 0.5 * (2 * y2 - yc1[ipad][irow][0] - yc2[ipad][irow][0]) * (x2 - x1);
            continue;
         }
         // case b) 0-1100
         if (indexPad[ipad][irow][0] == 0 && indexPad[ipad][irow][1] == 1100) {
            if (xc1[ipad][irow][1] < x2)
               sw1 = 0.5 * (x2 - xc1[ipad][irow][1]) * (y2 - yc2[ipad][irow][1]);
            else
               sw1 = 0.5 * (x2 - xc2[ipad][irow][1]) * (y2 - yc1[ipad][irow][1]);
            SPads[ipad] = SPad - sw1;
            continue;
         }
         // case d) 1100-0
         if (indexPad[ipad][irow][1] == 0 && indexPad[ipad][irow][0] == 1100) {
            if (xc1[ipad][irow][0] < x2)
               SPads[ipad] = 0.5 * (x2 - xc1[ipad][irow][0]) * (y2 - yc2[ipad][irow][0]);
            else
               SPads[ipad] = 0.5 * (x2 - xc2[ipad][irow][0]) * (y2 - yc1[ipad][irow][0]);
            continue;
         }
         //  * * * * * * * *
         // case f) 0011-0101
         if (indexPad[ipad][irow][0] == 11 && indexPad[ipad][irow][1] == 101) {
            if (yc1[ipad][irow][0] > y1)
               sw1 = 0.5 * (xc2[ipad][irow][0] - x1) * (yc1[ipad][irow][0] - y1);
            else
               sw1 = 0.5 * (xc1[ipad][irow][0] - x1) * (yc2[ipad][irow][0] - y1);
            sw2         = 0.5 * (2 * x2 - xc1[ipad][irow][1] - xc2[ipad][irow][1]) * hPad;
            SPads[ipad] = SPad - sw1 - sw2;
            continue;
         }
         // case f) 0101-0011
         if (indexPad[ipad][irow][1] == 11 && indexPad[ipad][irow][0] == 101) {
            if (yc1[ipad][irow][1] > y1)
               sw1 = 0.5 * (xc2[ipad][irow][1] - x1) * (yc1[ipad][irow][1] - y1);
            else
               sw1 = 0.5 * (xc1[ipad][irow][1] - x1) * (yc2[ipad][irow][1] - y1);
            sw2         = 0.5 * (2 * x2 - xc1[ipad][irow][0] - xc2[ipad][irow][0]) * hPad;
            SPads[ipad] = SPad - sw1 - sw2;
            continue;
         }
         // case g) 0011-0011
         if (indexPad[ipad][irow][0] == 11 && indexPad[ipad][irow][1] == 11) {
            if (yc1[ipad][irow][0] > y1)
               sw1 = 0.5 * (xc2[ipad][irow][0] - x1) * (yc1[ipad][irow][0] - y1);
            else
               sw1 = 0.5 * (xc1[ipad][irow][0] - x1) * (yc2[ipad][irow][0] - y1);
            if (yc1[ipad][irow][1] > y1)
               sw2 = 0.5 * (xc2[ipad][irow][1] - x1) * (yc1[ipad][irow][1] - y1);
            else
               sw2 = 0.5 * (xc1[ipad][irow][1] - x1) * (yc2[ipad][irow][1] - y1);
            SPads[ipad] = sw2 - sw1;
            continue;
         }
         // case g) 0011-1010
         if (indexPad[ipad][irow][0] == 11 && indexPad[ipad][irow][1] == 1010) {
            if (yc1[ipad][irow][0] > y1)
               sw1 = 0.5 * (xc2[ipad][irow][0] - x1) * (yc1[ipad][irow][0] - y1);
            else
               sw1 = 0.5 * (xc1[ipad][irow][0] - x1) * (yc2[ipad][irow][0] - y1);
            sw2         = 0.5 * (2 * y2 - yc1[ipad][irow][1] - yc2[ipad][irow][1]) * wPad;
            SPads[ipad] = SPad - sw1 - sw2;
            continue;
         }
         // case g) 0011-1100
         if (indexPad[ipad][irow][0] == 11 && indexPad[ipad][irow][1] == 1100) {
            if (yc1[ipad][irow][0] > y1)
               sw1 = 0.5 * (xc2[ipad][irow][0] - x1) * (yc1[ipad][irow][0] - y1);
            else
               sw1 = 0.5 * (xc1[ipad][irow][0] - x1) * (yc2[ipad][irow][0] - y1);
            if (yc1[ipad][irow][1] > y1)
               sw2 = 0.5 * (x2 - xc1[ipad][irow][1]) * (y2 - yc2[ipad][irow][1]);
            else
               sw2 = 0.5 * (x2 - xc2[ipad][irow][1]) * (y2 - yc1[ipad][irow][1]);
            SPads[ipad] = SPad - sw1 - sw2;
            continue;
         }
         // case g) 0011-1100
         if (indexPad[ipad][irow][0] == 1100 && indexPad[ipad][irow][1] == 11) {
            if (yc1[ipad][irow][1] > y1)
               sw1 = 0.5 * (xc2[ipad][irow][1] - x1) * (yc1[ipad][irow][1] - y1);
            else
               sw1 = 0.5 * (xc1[ipad][irow][1] - x1) * (yc2[ipad][irow][1] - y1);
            if (yc1[ipad][irow][0] > y1)
               sw2 = 0.5 * (x2 - xc1[ipad][irow][0]) * (y2 - yc2[ipad][irow][0]);
            else
               sw2 = 0.5 * (x2 - xc2[ipad][irow][0]) * (y2 - yc1[ipad][irow][0]);
            SPads[ipad] = SPad - sw1 - sw2;
            continue;
         }
         //  * * * * * * * *
         // case f) 0101-1001
         if (indexPad[ipad][irow][0] == 101 && indexPad[ipad][irow][1] == 1001) {
            if (yc1[ipad][irow][0] > y1)
               sw1 = 0.5 * (x2 - xc2[ipad][irow][0]) * (yc1[ipad][irow][0] - y1);
            else
               sw1 = 0.5 * (x2 - xc1[ipad][irow][0]) * (yc2[ipad][irow][0] - y1);
            sw2         = 0.5 * (xc1[ipad][irow][1] + xc2[ipad][irow][1] - 2 * x1) * hPad;
            SPads[ipad] = SPad - sw1 - sw2;
            continue;
         }
         // case f) 0101-0101
         if (indexPad[ipad][irow][0] == 101 && indexPad[ipad][irow][1] == 101) {
            SPads[ipad] = TMath::Abs(xc2[ipad][irow][0] - xc1[ipad][irow][1]) * hPad;
            continue;
         }
         // case f) 0101-0110
         if (indexPad[ipad][irow][0] == 101 && indexPad[ipad][irow][1] == 110) {
            if (xc1[ipad][irow][1] > x1)
               sw1 = 0.5 * (xc1[ipad][irow][1] - x1) * (y2 - yc2[ipad][irow][1]);
            else
               sw1 = 0.5 * (xc2[ipad][irow][1] - x1) * (y2 - yc1[ipad][irow][1]);
            sw2         = 0.5 * (yc1[ipad][irow][0] + yc2[ipad][irow][0] - 2 * y1) * wPad;
            SPads[ipad] = SPad - sw1 - sw2;
            continue;
         }
         // case f) 0110-0101
         if (indexPad[ipad][irow][1] == 101 && indexPad[ipad][irow][0] == 110) {
            if (xc1[ipad][irow][0] > x1)
               sw1 = 0.5 * (xc1[ipad][irow][0] - x1) * (y2 - yc2[ipad][irow][0]);
            else
               sw1 = 0.5 * (xc2[ipad][irow][0] - x1) * (y2 - yc1[ipad][irow][0]);
            sw2         = 0.5 * (yc1[ipad][irow][1] + yc2[ipad][irow][1] - 2 * y1) * wPad;
            SPads[ipad] = SPad - sw1 - sw2;
            continue;
         }
         // case f) 0101-1100
         if (indexPad[ipad][irow][0] == 101 && indexPad[ipad][irow][1] == 1100) {
            if (xc1[ipad][irow][1] < x2)
               sw1 = 0.5 * (x2 - xc1[ipad][irow][1]) * (y2 - yc2[ipad][irow][1]);
            else
               sw1 = 0.5 * (x2 - xc2[ipad][irow][1]) * (y2 - yc1[ipad][irow][1]);
            sw2         = 0.5 * (yc1[ipad][irow][0] + yc2[ipad][irow][0] - 2 * y1) * wPad;
            SPads[ipad] = SPad - sw1 - sw2;
            continue;
         }
         // case f) 1100-0101
         if (indexPad[ipad][irow][1] == 101 && indexPad[ipad][irow][0] == 1100) {
            if (xc1[ipad][irow][0] < x2)
               sw1 = 0.5 * (x2 - xc1[ipad][irow][0]) * (y2 - yc2[ipad][irow][0]);
            else
               sw1 = 0.5 * (x2 - xc2[ipad][irow][0]) * (y2 - yc1[ipad][irow][0]);
            sw2         = 0.5 * (yc1[ipad][irow][1] + yc2[ipad][irow][1] - 2 * y1) * wPad;
            SPads[ipad] = SPad - sw1 - sw2;
            continue;
         }
         //  * * * * * * * *
         // case g) 0101-1001
         if (indexPad[ipad][irow][0] == 1001 && indexPad[ipad][irow][1] == 1001) {
            if (yc1[ipad][irow][0] > y1)
               sw1 = 0.5 * (x2 - xc2[ipad][irow][0]) * (yc1[ipad][irow][0] - y1);
            else
               sw1 = 0.5 * (x2 - xc1[ipad][irow][0]) * (yc2[ipad][irow][0] - y1);
            if (yc1[ipad][irow][1] > y1)
               sw2 = 0.5 * (x2 - xc2[ipad][irow][1]) * (yc1[ipad][irow][1] - y1);
            else
               sw2 = 0.5 * (x2 - xc1[ipad][irow][1]) * (yc2[ipad][irow][1] - y1);
            SPads[ipad] = sw2 - sw1;
            continue;
         }
         // case f) 0101-0110
         if (indexPad[ipad][irow][0] == 101 && indexPad[ipad][irow][1] == 110) {
            if (xc1[ipad][irow][1] > x1)
               sw2 = 0.5 * (xc1[ipad][irow][1] - x1) * (y2 - yc2[ipad][irow][1]);
            else
               sw2 = 0.5 * (xc2[ipad][irow][1] - x1) * (y2 - yc1[ipad][irow][1]);
            sw1         = 0.5 * (2 * x2 - xc1[ipad][irow][0] - xc2[ipad][irow][0]) * hPad;
            SPads[ipad] = SPad - sw1 - sw2;
            continue;
         }
         // case e) 1001-0110
         if (indexPad[ipad][irow][0] == 1001 && indexPad[ipad][irow][1] == 110) {
            if (xc1[ipad][irow][0] < x2)
               sw1 = 0.5 * (x2 - xc1[ipad][irow][0]) * (yc2[ipad][irow][0] - y1);
            else
               sw1 = 0.5 * (x2 - xc2[ipad][irow][0]) * (yc1[ipad][irow][0] - y1);
            if (yc1[ipad][irow][1] < y2)
               sw2 = 0.5 * (xc2[ipad][irow][1] - x1) * (y2 - yc1[ipad][irow][1]);
            else
               sw2 = 0.5 * (xc1[ipad][irow][1] - x1) * (y2 - yc2[ipad][irow][1]);
            SPads[ipad] = SPad - sw1 - sw2;
            /*printf("r0c1(%f,%f)  r0c2(%f,%f)\n",
            xc1[ipad][irow][0],yc1[ipad][irow][0],xc2[ipad][irow][0],yc2[ipad][irow][0]);
            printf("r1c1(%f,%f)  r1c2(%f,%f)\n",
            xc1[ipad][irow][1],yc1[ipad][irow][1],xc2[ipad][irow][1],yc2[ipad][irow][1]);
            printf("SPads[%d]=%f-%f-%f=%f\n",ipad,SPad,sw1,sw2,SPads[ipad]);*/
            continue;
         }
         // case e) 110-1001
         if (indexPad[ipad][irow][0] == 110 && indexPad[ipad][irow][1] == 1001) {
            if (xc1[ipad][irow][1] < x2)
               sw1 = 0.5 * (x2 - xc1[ipad][irow][1]) * (yc2[ipad][irow][1] - y1);
            else
               sw1 = 0.5 * (x2 - xc2[ipad][irow][1]) * (yc1[ipad][irow][1] - y1);
            if (yc1[ipad][irow][0] < y2)
               sw2 = 0.5 * (xc2[ipad][irow][0] - x1) * (y2 - yc1[ipad][irow][0]);
            else
               sw2 = 0.5 * (xc1[ipad][irow][0] - x1) * (y2 - yc2[ipad][irow][0]);
            SPads[ipad] = SPad - sw1 - sw2;
            continue;
         }
         // case e) 1001-0101
         if (indexPad[ipad][irow][0] == 1001 && indexPad[ipad][irow][1] == 101) {
            if (xc1[ipad][irow][0] < x2)
               sw1 = 0.5 * (x2 - xc1[ipad][irow][0]) * (yc2[ipad][irow][0] - y1);
            else
               sw1 = 0.5 * (x2 - xc2[ipad][irow][0]) * (yc1[ipad][irow][0] - y1);
            sw2         = 0.5 * (2 * y2 - yc1[ipad][irow][1] - yc2[ipad][irow][1]) * wPad;
            SPads[ipad] = SPad - sw1 - sw2;
            continue;
         }
         // case e) 1001-1010
         if (indexPad[ipad][irow][0] == 1001 && indexPad[ipad][irow][1] == 1010) {
            if (xc1[ipad][irow][0] < x2)
               sw1 = 0.5 * (x2 - xc1[ipad][irow][0]) * (yc2[ipad][irow][0] - y1);
            else
               sw1 = 0.5 * (x2 - xc2[ipad][irow][0]) * (yc1[ipad][irow][0] - y1);
            sw2         = 0.5 * (2 * y2 - yc1[ipad][irow][1] - yc2[ipad][irow][1]) * wPad;
            SPads[ipad] = SPad - sw1 - sw2;
            continue;
         }
         // case g) 1010-1010
         if (indexPad[ipad][irow][0] == 1010 && indexPad[ipad][irow][1] == 1010) {
            SPads[ipad] =
               0.5 * (yc1[ipad][irow][1] - yc1[ipad][irow][0] + yc2[ipad][irow][1] - yc2[ipad][irow][0]) * wPad;
            continue;
         }
         // case g) 1010-0110
         if (indexPad[ipad][irow][0] == 1010 && indexPad[ipad][irow][1] == 110) {
            if (xc1[ipad][irow][1] > x1)
               sw1 = 0.5 * (xc1[ipad][irow][1] - x1) * (y2 - yc2[ipad][irow][1]);
            else
               sw1 = 0.5 * (xc2[ipad][irow][1] - x1) * (y2 - yc1[ipad][irow][1]);
            sw2         = 0.5 * (yc1[ipad][irow][0] + yc2[ipad][irow][0] - 2 * y1) * wPad;
            SPads[ipad] = SPad - sw1 - sw2;
            /*printf("1010-0110 r0c1(%f,%f)  r0c2(%f,%f)\n",
            xc1[ipad][irow][0],yc1[ipad][irow][0],xc2[ipad][irow][0],yc2[ipad][irow][0]);
            printf("1010-0110 r1c1(%f,%f)  r1c2(%f,%f)\n",
            xc1[ipad][irow][1],yc1[ipad][irow][1],xc2[ipad][irow][1],yc2[ipad][irow][1]);
            printf("SPads[%d]=%f-%f-%f=%f\n",ipad,SPad,sw1,sw2,SPads[ipad]);*/
            continue;
         }
         // case g) 1010-1100
         if (indexPad[ipad][irow][0] == 1010 && indexPad[ipad][irow][1] == 1100) {
            if (xc1[ipad][irow][1] < x2)
               sw1 = 0.5 * (x2 - xc1[ipad][irow][1]) * (y2 - yc2[ipad][irow][1]);
            else
               sw1 = 0.5 * (x2 - xc2[ipad][irow][1]) * (y2 - yc1[ipad][irow][1]);
            sw2         = 0.5 * (yc1[ipad][irow][0] + yc2[ipad][irow][0] - 2 * y1) * wPad;
            SPads[ipad] = SPad - sw1 - sw2;
            continue;
         }
         // case g) 0110-0110
         if (indexPad[ipad][irow][0] == 110 && indexPad[ipad][irow][1] == 110) {
            if (yc1[ipad][irow][0] < y2)
               sw1 = 0.5 * (xc2[ipad][irow][0] - x1) * (y2 - yc1[ipad][irow][0]);
            else
               sw1 = 0.5 * (xc1[ipad][irow][0] - x1) * (y2 - yc2[ipad][irow][0]);
            if (yc1[ipad][irow][1] < y2)
               sw2 = 0.5 * (xc2[ipad][irow][1] - x1) * (y2 - yc1[ipad][irow][1]);
            else
               sw2 = 0.5 * (xc1[ipad][irow][1] - x1) * (y2 - yc2[ipad][irow][1]);
            SPads[ipad] = sw1 - sw2;
            continue;
         }
         // case g) 1100-1100
         if (indexPad[ipad][irow][0] == 1100 && indexPad[ipad][irow][1] == 1100) {
            if (yc1[ipad][irow][0] < y2)
               sw1 = 0.5 * (x2 - xc2[ipad][irow][0]) * (y2 - yc1[ipad][irow][0]);
            else
               sw1 = 0.5 * (x2 - xc1[ipad][irow][0]) * (y2 - yc2[ipad][irow][0]);
            if (yc1[ipad][irow][1] < y2)
               sw2 = 0.5 * (x2 - xc2[ipad][irow][1]) * (y2 - yc1[ipad][irow][1]);
            else
               sw2 = 0.5 * (x2 - xc1[ipad][irow][1]) * (y2 - yc2[ipad][irow][1]);
            SPads[ipad] = sw1 - sw2;
            /*printf("1100-1100 r0c1(%f,%f)  r0c2(%f,%f)\n",
            xc1[ipad][irow][0],yc1[ipad][irow][0],xc2[ipad][irow][0],yc2[ipad][irow][0]);
            printf("          r1c1(%f,%f)  r1c2(%f,%f)\n",
            xc1[ipad][irow][1],yc1[ipad][irow][1],xc2[ipad][irow][1],yc2[ipad][irow][1]);
            printf("            r1(%f,%f)    r2(%f,%f)\n",x1,y1,x2,y2);
            printf("SPads[%d]=%f-%f-%f=%f\n",ipad,SPad,sw1,sw2,SPads[ipad]);*/
            continue;
         }
         SPads[ipad] = 0;

         printf("case g110 rp[%d,%d] index(%d %d %d)\n", irow, ipad, indexPad[ipad][irow][0], indexPad[ipad][irow][1],
                indexPad[ipad][irow][2]);
         printf("r1(%f,%f) r2(%f,%f) \nr0c1(%f,%f) r0c2(%f,%f)\n", x1, y1, x2, y2, xc1[ipad][irow][0],
                yc1[ipad][irow][0], xc2[ipad][irow][0], yc2[ipad][irow][0]);
         printf("r1c1(%f,%f) r1c2(%f,%f)\n", xc1[ipad][irow][1], yc1[ipad][irow][1], xc2[ipad][irow][1],
                yc2[ipad][irow][1]);
         // *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *

      } // end of pads cycle
        /*for(Int_t ip=ip1;ip<=ip2;ip++) {
        printf(" p[%d %d %d] tl0(%f-%f) tl1(%f-%f) tl2(%f-%f) indx(%d %d %d)\n",
        iSector,ip,irow,tPad1[ip][irow][0],
        tPad2[ip][irow][0],tPad1[ip][irow][1],tPad1[ip][irow][1],tPad1[ip][irow][2],
        tPad1[ip][irow][2],indexPad[ip][irow][0],indexPad[ip][irow][1],indexPad[ip][irow][2]);
        }*/
        //   ***************************************************
        //   ********       check pads in each row      ********
        //   ***************************************************
      Int_t iRoc = fTpcSecGeo->NofRoc(irow);
      //        Int_t iP2=ip2-1;
      nclust += ClustersInRow(irow, ip1, ip2, iRoc, SPads);
   } // end of row  cycle
}

//______________________________________________________________________
// Get get helix point coordinates and sector number
Int_t TpcMissAlignment::coordinates(Double_t t, Double_t R)
{
   Double_t phi;
   Int_t    sec;
   rBt = TVector3(xcB + R * TMath::Cos(phi_0B + signRotB * t), // BCS
                  ycB + R * TMath::Sin(phi_0B + signRotB * t), r0BeamB.Z() + helix_Pz * t);
   // printf("rBt(%f,%f,%f) rcB(%f,%f) phi_0B=%f t=%f\n",rBt.X(),rBt.Y(),rBt.Z(),
   //         xcB,ycB,phi_0B*r2d,t*r2d);
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
      R0shift  = rR0_A[sec];
      loc2glo  = r_loc2glo[sec];
      e_3      = TVector3(loc2glo.XZ(), loc2glo.YZ(), loc2glo.ZZ());
      G0shift  = rG0_A[sec];
      glo2loc  = r_glo2loc[sec];
      uR0shift = uR0_A[sec];
      uloc2glo = u_loc2glo[sec];
      e_3u     = TVector3(u_loc2glo[sec].XZ(), u_loc2glo[sec].YZ(), u_loc2glo[sec].ZZ());
   }
   // printf("rBt(%f,%f,%f) rcB(%f,%f) phi_0B=%f  sec=%d\n",rBt.X(),rBt.Y(),rBt.Z(),
   //         xcB,ycB,phi_0B*r2d,sec);
   Double_t t_s   = (e_3 * (R0shift - rGt)) / (e_3 * eEf); // t to sector plane point
   TVector3 rGt_s = rGt + t_s * eEf;                       // helix point on the sector plane (GCS)
   rLt            = G0shift + glo2loc * rGt_s;
   // printf("rLt(%f,%f,%f) \n",rLt.X(),rLt.Y(),rLt.Z());
   return sec;
}

//______________________________________________________________________
// Get get helix point coordinates and sector number
void TpcMissAlignment::coordinates(Int_t sec, Double_t t, Double_t R)
{
   Double_t phi;
   rBt = TVector3(xcB + R * TMath::Cos(phi_0B + signRotB * t), // BCS
                  ycB + R * TMath::Sin(phi_0B + signRotB * t), r0BeamB.Z() + helix_Pz * t);
   // printf("rBt(%f,%f,%f) rcB(%f,%f) phi_0B=%f\n",rBt.X(),rBt.Y(),rBt.Z(),
   //         xcB,ycB,phi_0B*r2d);
   rGt = B2G * rBt;
   // printf("rBt(%f,%f,%f) rcB(%f,%f) phi_0B=%f  sec=%d\n",rBt.X(),rBt.Y(),rBt.Z(),
   //         xcB,ycB,phi_0B*r2d,sec);
   Double_t t_s   = (e_3 * (R0shift - rGt)) / (e_3 * eEf); // t to sector plane point
   TVector3 rGt_s = rGt + t_s * eEf;                       // helix point on the sector plane (GCS)
   rLt            = G0shift + glo2loc * rGt_s;
   // printf("rLt(%f,%f,%f) \n",rLt.X(),rLt.Y(),rLt.Z());
   return;
}

//______________________________________________________________________

// simulate z pozition
Double_t TpcMissAlignment::simHelixZpos(TVector2 LCxy, Int_t ipad)
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
   Double_t t = tPads[ipad];
   // find particle coordinates in BCS
   TVector3 particleB = TVector3(xcB + rH * TMath::Cos(phi_0B + signRotB * t), // BCS
                                 ycB + rH * TMath::Sin(phi_0B + signRotB * t), r0BeamB.Z() + helix_Pz * t);
   // find particle coordinates in GCS
   TVector3 particleG = B2G * particleB;
   // find the distance (D) from the particle to the sector plane
   TVector3 D = SecB - particleB; // vector: cluster->track point
   Double_t z = fRandom->Gaus(D.Mag(), D.Mag() * sig_Z_pos);
   if (z < 0) {
      z = 0;
   } else {
      if (z > fL_svol) {
         z = sZ_min;
      }
   }
   // printf("simHelixZpos: t=%f rL(%f,%f) rpG(%f,%f,%f)\n",t,LCxy.X(),LCxy.Y(),
   // particleG.X(),particleG.Y(),
   // particleG.Z());
   // printf("D.Mag()=%f z=%f sZ_min=%f\n",D.Mag(),z,sZ_min);
   // printf("simHelixZpos:t=%f particleB(%f,%f,%f) d=%f z=%f\n",
   // t,particleB.X(),particleB.Y(),particleB.Z(),D.Mag(),z);
   return z;
}

//______________________________________________________________________
// Chi2 for the circle projection
Double_t TpcMissAlignment_Chi2Circle(const Double_t *par)
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
   for (int i = 0; i < MpdTpcMiniMC_Nhits; i++) {
      X          = MpdTpcMiniMC_XBHits[i];
      Y          = MpdTpcMiniMC_YBHits[i];
      double dfi = phi0 - TMath::ATan2(Y - yc, X - xc);
      if (MpdTpcMiniMC_qHelix > 0) {
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
      double dx = X - (xc + R * TMath::Cos(phi0 - MpdTpcMiniMC_qHelix * t));
      double dy = Y - (yc + R * TMath::Sin(phi0 - MpdTpcMiniMC_qHelix * t));
      // printf("i=%2d t=%6.1f tF=%6.1f H(%6.1f %6.1f) T(%6.1f %6.1f) dx=%f dy=%f chi2=%f\n",i,
      // t*57.29577951,(phi0-MpdTpcMiniMC_qHelix*t)*57.29577951,X,Y,X-dx,Y-dy,dx,dy,chisq);
      delta = TMath::Sqrt(dx * dx + dy * dy) / MpdTpcMiniMC_WHits[i];
      chisq += delta * delta;
   }
   // printf("xc,yc(%6.1f %6.1f) R=%6.1f phi0=%6.1f chi2=%f\n",
   // xc,yc,R,phi0*57.29577951,chisq);
   return chisq;
}

//______________________________________________________________________
// Chi2 for the clusters fit by the helix
Double_t TpcMissAlignment_zHelix(const Double_t *par)
{

   double chisq = 0;
   double xc    = trackPout[0];
   double yc    = trackPout[1];
   double R     = TMath::Abs(trackPout[3]);
   double phi0  = trackPout[4];
   // printf("xc,yc(%6.1f %6.1f) R=%6.1f phi0=%6.1f\n",xc,yc,R,phi0*57.29577951);
   for (int i = 0; i < MpdTpcMiniMC_Nhits; i++) {
      double t     = MpdTpcMiniMC_tHelix[i];
      double dx    = MpdTpcMiniMC_XBHits[i] - (xc + R * TMath::Cos(phi0 - MpdTpcMiniMC_qHelix * t));
      double dy    = MpdTpcMiniMC_YBHits[i] - (yc + R * TMath::Sin(phi0 - MpdTpcMiniMC_qHelix * t));
      double dz    = (MpdTpcMiniMC_ZBHits[i] - (par[0] + par[1] * MpdTpcMiniMC_tHelix[i]));
      double delta = TMath::Sqrt(dx * dx + dy * dy + dz * dz) / MpdTpcMiniMC_WHits[i];
      // printf("i=%2d t=%6.1f (%8.1f %8.1f) dx=%f dy=%f dz=%f chi2=%f\n",i,
      // MpdTpcMiniMC_tHelix[i],par[0],par[1],dx,dy,dz,chisq);
      chisq += delta * delta;
   }
   // printf("par(%lf %lf) chi2=%lf\n",par[0],par[1],chisq);
   return chisq;
}

//______________________________________________________________________
// Chi2 for the clusters fit by the helix
Double_t TpcMissAlignment_Chi2Helix(const Double_t *par)
{
   double t;
   double xc     = par[0];
   double yc     = par[1];
   double zc     = par[2];
   double R      = par[3];
   double phi0   = par[4];
   double zHelix = par[5];
   int    q      = MpdTpcMiniMC_qHelix;
   double X      = xc + R * TMath::Cos(phi0);
   double Y      = yc + R * TMath::Sin(phi0);
   double delta  = TMath::Sqrt(X * X + Y * Y) / 0.2;
   double chisq  = delta * delta;
   for (int i = 0; i < MpdTpcMiniMC_Nhits; i++) {
      X          = MpdTpcMiniMC_XBHits[i];
      Y          = MpdTpcMiniMC_YBHits[i];
      double Z   = MpdTpcMiniMC_ZBHits[i];
      double dfi = phi0 - TMath::ATan2(Y - yc, X - xc);
      if (q > 0) {
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
      double dx = X - (xc + R * TMath::Cos(phi0 - q * t));
      double dy = Y - (yc + R * TMath::Sin(phi0 - q * t));
      double dz = Z - (zc + zHelix * t);
      delta     = TMath::Sqrt(dx * dx + dy * dy + dz * dz) / MpdTpcMiniMC_WHits[i];
      chisq += delta * delta;
   }
   return chisq;
}
// printf("   chi2[%d]=%f\n",i,chisq);

//______________________________________________________________________
// Reconstruction of the helix parameters
Int_t TpcMissAlignment::fitHelix()
{
   /*
     in the BCS the helix equation:
       x=p[0]+p[3]*cos(p[4]+q*t)
       y=p[1]+p[3]*sin(p[4]+q*t)
       z=p[2]+p[5]*t
   */
   Int_t    imin, imax;
   Double_t zmin = 200, zmax = 0;
   for (Int_t i = 0; i < MpdTpcMiniMC_Nhits; i++) {
      if (MpdTpcMiniMC_ZBHits[i] < zmin) {
         imin = i;
         zmin = MpdTpcMiniMC_ZBHits[i];
      }
      if (MpdTpcMiniMC_ZBHits[i] > zmax) {
         imax = i;
         zmax = MpdTpcMiniMC_ZBHits[i];
      }
      // printf("i=%2d sec=%2d XYZ(%8.2f %8.2f %8.2f)\n",i,MpdTpcMiniMC_iSec[i],
      // MpdTpcMiniMC_XBHits[i],MpdTpcMiniMC_YBHits[i],MpdTpcMiniMC_ZBHits[i]);
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
   x2 = MpdTpcMiniMC_XBHits[imin];
   y2 = MpdTpcMiniMC_YBHits[imin];
   z2 = MpdTpcMiniMC_ZBHits[imin];
   x3 = MpdTpcMiniMC_XBHits[imax];
   y3 = MpdTpcMiniMC_YBHits[imax];
   z3 = MpdTpcMiniMC_ZBHits[imax];

   CircleCenter(x1, y1, x2, y2, x3, y3, xc, yc),
      //    xc+=xw; yc+=yw;
      //    CircleCenter(x2,y2,x3,y3,x1,y1,xw,yw),
      //    xc+=xw; yc+=yw;
      //    CircleCenter(x3,y3,x1,y1,x2,y2,xw,yw),
      //    xc+=xw; yc+=yw;
      //    xc/=3; yc/=3;*/

      // determ the particle sign
      t1 = TMath::ATan2(y2 - yc, x2 - xc);
   t2    = TMath::ATan2(y3 - yc, x3 - xc);
   if (t1 > t2) {
      if (t1 - t2 > pi) {
         MpdTpcMiniMC_qHelix = -1;
         t                   = twopi - (t1 - t2);
      } else {
         MpdTpcMiniMC_qHelix = 1;
         t                   = t1 - t2;
      }
   } else {
      if (t2 - t1 > pi) {
         MpdTpcMiniMC_qHelix = 1;
         t                   = twopi - (t2 - t1);
      } else {
         MpdTpcMiniMC_qHelix = -1;
         t                   = t2 - t1;
      }
   }
   // printf("fitHelix: R1(%8.1f,%8.1f) R2(%8.1f %8.1f)  q=%d\n",x2,y2,x3,y3,MpdTpcMiniMC_qHelix);

   // Choose method upon creation between:
   // kMigrad, kSimplex, kCombined,
   // kScan, kFumili
   static ROOT::Minuit2::Minuit2Minimizer minFitc(ROOT::Minuit2::kMigrad);
   //  ================ find 4 circle parameters =================
   ROOT::Math::Functor fcn1(&TpcMissAlignment_Chi2Circle, 4);
   //                     Xc    Yc    R     Phi0
   double step[6] = {1, 1, 1, 0.02, 0, 0};
   double limU[6] = {2000, 2000, 2000, TMath::Pi(), 0, 0};
   double limL[6] = {-2000, -2000, -2000, -TMath::Pi(), 0, 0};
   double par[6];

   par[0] = xc;                             //
   par[1] = yc;                             //  - coordinates of the muon start point
   par[2] = TMath::Sqrt(xc * xc + yc * yc); // radius of the circle in XY plane of BCS
   par[3] = TMath::ATan2(-yc, -xc);         // initial muon angle (XY plane of BCS)

   /*printf("PinHelix: par(%8.1f %8.1f %8.1f %6.1f)  q=%d nhits=%d\n",
   trackPin[0],trackPin[1],trackPin[3],r2d*trackPin[4],
   MpdTpcMiniMC_qHelix,MpdTpcMiniMC_Nhits);
   printf("fitHelix: par(%8.1f %8.1f %8.1f %6.1f)  q=%d\n",
   par[0],par[1],par[2],r2d*par[3],MpdTpcMiniMC_qHelix);*/

   minFitc.SetFunction(fcn1);

   // Set the free variables to be minimized!
   minFitc.SetLimitedVariable(0, "x0", par[0], step[0], limL[0], limU[0]);
   minFitc.SetLimitedVariable(1, "y0", par[1], step[1], limL[1], limU[1]);
   minFitc.SetLimitedVariable(2, "R", par[2], step[2], limL[2], limU[2]);
   minFitc.SetLimitedVariable(3, "phi", par[3], step[3], limL[3], limU[3]);

   Bool_t        res = minFitc.Minimize();
   const double *xs1 = minFitc.X();
   /*printf("outHelix: par(%8.1f %8.1f %8.1f %6.1f)  q=%d\n",
   xs1[0],xs1[1],xs1[2],r2d*xs1[3],MpdTpcMiniMC_qHelix);
   printf("chi2min4=%f\n",minFitc.MinValue()/MpdTpcMiniMC_Nhits);
   cout<<"res="<<res<<endl<<endl;*/
   if (minFitc.MinValue() / MpdTpcMiniMC_Nhits > EvalpHit) return 1;
   //      minFitc.PrintResults();
   xc = trackPout[0] = xs1[0];
   yc = trackPout[1] = xs1[1];
   trackPout[3]      = MpdTpcMiniMC_qHelix * xs1[2];
   Double_t phi0 = trackPout[4] = xs1[3];

   trackDiff[0] = trackPout[0] - trackPin[0];
   trackDiff[1] = trackPout[1] - trackPin[1];
   trackDiff[3] = trackPout[3] - trackPin[3];
   trackDiff[4] = trackPout[4] - trackPin[4];

   //  ================ find 2 longitudal parameters ================
   static ROOT::Minuit2::Minuit2Minimizer minFitz(ROOT::Minuit2::kMigrad);

   Double_t R = TMath::Abs(trackPout[3]);
   for (int i = 0; i < MpdTpcMiniMC_Nhits; i++) {
      Double_t X   = MpdTpcMiniMC_XBHits[i];
      Double_t Y   = MpdTpcMiniMC_YBHits[i];
      Double_t phi = TMath::ATan2(Y - yc, X - xc);
      if (MpdTpcMiniMC_qHelix > 0)
         MpdTpcMiniMC_tHelix[i] =
            (phi > phi0) ? twopi - MpdTpcMiniMC_qHelix * (phi - phi0) : -MpdTpcMiniMC_qHelix * (phi - phi0);
      else
         MpdTpcMiniMC_tHelix[i] =
            (phi > phi0) ? -MpdTpcMiniMC_qHelix * (phi - phi0) : twopi - MpdTpcMiniMC_qHelix * (phi - phi0);
      /*printf("i=%2d R=(%6.1f %6.1f %6.1f) t=%6.1f phi0=%7.1f phi=%7.1f\n",i,
      MpdTpcMiniMC_XBHits[i],MpdTpcMiniMC_YBHits[i],MpdTpcMiniMC_ZBHits[i],
      MpdTpcMiniMC_tHelix[i],phi0*r2d,phi*r2d);*/
   }
   ROOT::Math::Functor fcn2(&TpcMissAlignment_zHelix, 2);
   step[0] = 0.1;
   step[1] = 0.1;
   limU[0] = 100;
   limU[1] = 2000;
   limL[0] = -100;
   limL[1] = -2000;

   par[0] = 0; //
   par[1] = (MpdTpcMiniMC_ZBHits[imax] - MpdTpcMiniMC_ZBHits[imin]) /
            (MpdTpcMiniMC_tHelix[imax] - MpdTpcMiniMC_tHelix[imin]); // initial Z-velocity
   /*printf("t1,t2,t(%6.1f %6.1f %6.1f) z2,z3(%6.1f %6.1f)\n",t1,t2,t,z2,z3);
   printf("PinHelix: par(%8.1f %8.1f)\n",trackPin[2],trackPin[5]);
   printf("fitHelix: par(%8.1f %8.1f)\n",par[0],par[1]);*/

   minFitz.SetFunction(fcn2);

   // Set the free variables to be minimized!
   minFitz.SetLimitedVariable(0, "x0", par[0], step[0], limL[0], limU[0]);
   minFitz.SetLimitedVariable(1, "y0", par[1], step[1], limL[1], limU[1]);

   res               = minFitz.Minimize();
   const double *xs2 = minFitz.X();
   if (minFitz.MinValue() / MpdTpcMiniMC_Nhits > EvalpHit) return 2;
   /*printf("outHelix: par(%8.1f %8.1f)\n",xs2[0],xs2[1]);
   printf("chi2min2=%f\n",minFitz.MinValue()/MpdTpcMiniMC_Nhits);
         minFitz.PrintResults();
   cout<<"res="<<res<<endl;*/
   trackPout[2] = xs2[0];
   trackPout[5] = xs2[1];

   /*
   //  ================ find 6 longitudal parameters ================
       static ROOT::Minuit2::Minuit2Minimizer minFit(ROOT::Minuit2::kMigrad);
       ROOT::Math::Functor fcn(&TpcMissAlignment_Chi2Helix,6);

       step[0]=1;step[1]=0.1;step[2]=0.1;step[3]=0.1;step[4]=0.02;step[5]=0.1;
       limU[0]=2000;limU[1]=2000;limU[2]=2000;limU[3]=2000;
       limU[4]=TMath::Pi();limU[5]=2000;
       limL[0]=-2000;limL[1]=-2000;limL[2]=-2000;limL[3]=-2000;
       limL[4]=-TMath::Pi();limL[5]=-2000;

       for(Int_t i=0;i<6;i++) par[i]=trackPout[i];

   printf("PinHelix: par(%8.1f %8.1f %4.1f %8.1f %6.1f %8.1f) q=%d nhits=%d\n",
   trackPin[0],trackPin[1],trackPin[2],trackPin[3],r2d*trackPin[4],trackPin[5],
   MpdTpcMiniMC_qHelix,MpdTpcMiniMC_Nhits);
   printf("fitHelix: par(%8.1f %8.1f %4.1f %8.1f %6.1f %8.1f)  q=%d\n\n",
   par[0],par[1],par[2],par[3],r2d*par[4],par[5],MpdTpcMiniMC_qHelix);

       minFit.SetFunction(fcn);

       // Set the free variables to be minimized!
       minFit.SetLimitedVariable(0,"x0",   par[0],step[0],limL[0],limU[0]);
       minFit.SetLimitedVariable(1,"y0",   par[1],step[1],limL[1],limU[1]);
       minFit.SetLimitedVariable(2,"z0",   par[2],step[2],limL[2],limU[2]);
       minFit.SetLimitedVariable(3,"R",    par[3],step[3],limL[3],limU[3]);
       minFit.SetLimitedVariable(4,"phi",  par[4],step[4],limL[4],limU[4]);
       minFit.SetLimitedVariable(5,"Pz",   par[5],step[5],limL[5],limU[5]);

       res=minFit.Minimize();
         const double *xs = minFit.X();
   printf("outHelix: par(%8.1f %8.1f %4.1f %8.1f %6.1f %8.1f)  q=%d\n",
   xs[0],xs[1],xs[2],xs[3],r2d*xs[4],xs[5],MpdTpcMiniMC_qHelix);
   printf("chi2min6=%f\n",minFit.MinValue()/MpdTpcMiniMC_Nhits);
   cout<<"res="<<res<<endl;
   //      if(minFit.MinValue()/MpdTpcMiniMC_Nhits>EvalpHit) return 1;
         minFit.PrintResults();
         trackPout[0]=xs[0];
         trackPout[1]=xs[1];
         trackPout[2]=xs[2];
         trackPout[3]=MpdTpcMiniMC_qHelix*xs[3];
         trackPout[4]=xs[4];
         trackPout[5]=xs[5];
   */

   for (Int_t k = 0; k < 6; k++) trackDiff[k] = trackPout[k] - trackPin[k];

   //      pRecovery();
   return 0;
}

//______________________________________________________________________
void TpcMissAlignment::CircleCenter(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Double_t x3, Double_t y3,
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
void TpcMissAlignment::putChi2H()
{
   //    for(Int_t i=0;i<6;i++) trackPout[i]=trackPin[i];
   Double_t chisq = 0;
   Double_t delta;
   Double_t phi0 = trackPout[4];
   Double_t R    = TMath::Abs(trackPout[3]);
   Int_t    q    = trackPout[3] / R;
   for (int i = 0; i < MpdTpcMiniMC_Nhits; i++) {
      Double_t t  = MpdTpcMiniMC_tHelix[i];
      Double_t dx = MpdTpcMiniMC_XBHits[i] - (trackPout[0] + R * TMath::Cos(phi0 - q * t));
      Double_t dy = MpdTpcMiniMC_YBHits[i] - (trackPout[1] + R * TMath::Sin(phi0 - q * t));
      Double_t dz = MpdTpcMiniMC_ZBHits[i] - (trackPout[2] + trackPout[5] * t);
      delta       = TMath::Sqrt(dx * dx + dy * dy + dz * dz) / MpdTpcMiniMC_WHits[i];
      chisq += delta * delta;
      MpdTpcMiniMC_chi2s[i] = delta * delta;
      fNhitSec[MpdTpcMiniMC_iSec[i]]++;
      fSecChi2[MpdTpcMiniMC_iSec[i]] += delta * delta;
   }
   Chi2Fit = chisq / MpdTpcMiniMC_Nhits;
}

//______________________________________________________________________
// momentum Recovery from helix parameters
void TpcMissAlignment::pRecovery()
{

   Double_t pT = c_p_m1 * TMath::Abs(trackPout[3]);
   printf("recovered pT=%f  c_p_m1=%f rH=%f\n", pT, c_p_m1, TMath::Abs(trackPout[3]));
   Double_t signR = (trackPout[3] > 0) ? -1 : 1;
   TVector3 pB    = TVector3(-signR * TMath::Sin(trackPout[4]) * pT, signR * TMath::Cos(trackPout[4]) * pT,
                             trackPout[5] * pT / TMath::Abs(trackPout[3]));
   TVector3 p     = B2G * pB;
   printf("input     p(%f,%f,%f)  charg=", pLas.X(), pLas.Y(), pLas.Z());
   if (q_sign > 0)
      printf("+1\n");
   else
      printf("-1\n");
   printf("recovered p(%f,%f,%f)  charg=", p.X(), p.Y(), p.Z());
   if (signR < 0)
      printf("+1\n");
   else
      printf("-1\n");
   return;
}

//______________________________________________________________________
// Get a ray of the TPC Laser System
void TpcMissAlignment::GetLaserRay()
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
void TpcMissAlignment::setMuon(Double_t x, Double_t y, Double_t z, Double_t theta, Double_t phi)
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
void TpcMissAlignment::GetBeamMuon()
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
// service subcode for pads in a row crossed by 2 circles: "/.../" or "\...\"
void TpcMissAlignment::n0is1n1is1(Int_t irow, Int_t ipl0[2][2], Int_t ipl1[2][2], Int_t &ip1, Int_t &ip2)
{
   Int_t    ipw1, ipw2;
   Double_t t1, t2;
   //  start between l0&l1
   // printf("n0is1n1is1: row=%d ipl0[%d-%d] ipl1[%d-%d]\n",irow,ipl0[0][0],ipl0[0][1],
   // ipl1[0][0],ipl1[0][1]);
   if (ipl0[0][0] <= ipl1[0][0]) {
      ip1 = ipl0[0][0];
      t1  = 0.5 * (tPad1[ipl0[0][1]][irow][0] + tPad2[ipl0[0][1]][irow][0]);
      t2  = 0.5 * (tPad1[ipl1[0][1]][irow][1] + tPad2[ipl1[0][1]][irow][1]);
      if (ipl1[0][1] >= ipl0[0][1])
         ip2 = ipl1[0][1]; //   I-l0-....l1-I
      else
         ip2 = ipl0[0][1];
   } else {
      ip1 = ipl1[0][0];
      t1  = 0.5 * (tPad1[ipl1[0][1]][irow][0] + tPad2[ipl1[0][1]][irow][0]);
      t2  = 0.5 * (tPad1[ipl0[0][1]][irow][1] + tPad2[ipl0[0][1]][irow][1]);
      if (ipl0[0][1] >= ipl1[0][1])
         ip2 = ipl0[0][1];
      else
         ip2 = ipl1[0][1];
   }
   // printf("n0is1n1is1: ip1=%d ip2=%d\n",ip1,ip2);

   for (Int_t ip = ip1; ip <= ip2; ip++) {
      if (indexPad[ip][irow][0] != 0) {
         tPads[ip] = 0.5 * (tPad1[ip][irow][0] + tPad2[ip][irow][0]);
      } else {
         if (indexPad[ip][irow][1] != 0) {
            tPads[ip] = 0.5 * (tPad1[ip][irow][1] + tPad2[ip][irow][1]);
         } else {
            indexPad[ip][irow][2] = 1111;
            tPads[ip]             = 0.5 * (t1 + t2);
         }
      }
      // printf("n0is1n1is1: indexPad[%d]=(%d %d %d) t=%f\n",ip,indexPad[ip][irow][0],
      // indexPad[ip][irow][1],indexPad[ip][irow][2],tPads[ip]);
   }
}

//______________________________________________________________________
// service subcode for pads in a row crossed by 1 circles: "/" or "\"

void TpcMissAlignment::n0orn1is1(Int_t n0, Int_t n1, Int_t irow, Int_t nPads, Int_t ipl0[2][2], Int_t ipl1[2][2],
                                 Int_t &ip1, Int_t &ip2)
{
   Int_t    ipw1, ipw2, ipt;
   Double_t t;
   if (n0 > 0) {
      ip1 = ipl0[0][0];
      ip2 = ipl0[0][1];
      // printf("n0 n0orn1is1: rp[%d,%d-%d] \n",irow,ip1,ip2);
      // printf("n0 n0orn1is1: ip1(0-1)(%f-%f)  ip2(0-1)(%f-%f)\n",tPad1[ip1][irow][0],
      // tPad2[ip1][irow][0],tPad1[ip2][irow][0],tPad2[ip2][irow][0]);
      if (xc1[ip2][irow][0] < xcL) { // signal pads on left
                                     // printf("n0orn1is1: n0left\n");
         ip1 = 0;
         if (yc1[ip1][irow][0] < yc2[ip1][irow][0]) // motion Up
            t = tPad1[ip1][irow][0];
         else
            t = tPad2[ip1][irow][0]; // motion Down
      } else {
         if (xc1[ip1][irow][0] > xcL) { // signal pads on right
                                        // printf("n0orn1is1: n0right\n");
            ip2 = nPads - 1;
            if (yc1[ip1][irow][0] > yc2[ip1][irow][0]) // motion Down
               t = tPad2[ip2][irow][0];
            else // motion Up
               t = tPad1[ip2][irow][0];
            // printf("n0orn1is1: n0rightUp\n");
         }
      }
   } else { // n1>0
      // printf("n1 n0orn1is1: rp[%d,%d-%d] \n",irow,ip1,ip2);
      ip1 = ipl1[0][0];
      ip2 = ipl1[0][1];
      if (xc1[ip2][irow][0] < xcL) { // signal pads on right
         ip2 = nPads - 1;
         if (yc1[ip1][irow][0] < yc2[ip1][irow][0]) // motion Up
            t = tPad2[ip2][irow][1];
         else // motion Down
            t = tPad1[ip2][irow][1];
      } else {
         if (xc1[ip1][irow][0] > xcL) { // signal pads on left
            ip1 = 0;
            if (yc1[ip1][irow][0] > yc2[ip1][irow][0]) // motion Down
               t = tPad1[ip1][irow][1];                // motion Down
            else
               t = tPad2[ip1][irow][1]; // motion Up
         }
      }
   }

   for (Int_t ip = ip1; ip <= ip2; ip++) {
      if (indexPad[ip][irow][0] != 0) {
         tPads[ip] = 0.5 * (tPad1[ip][irow][0] + tPad1[ip][irow][0]);
      } else {
         if (indexPad[ip][irow][1] != 0) {
            tPads[ip] = 0.5 * (tPad1[ip][irow][1] + tPad1[ip][irow][1]);
         } else {
            indexPad[ip][irow][2] = 1111;
            tPads[ip]             = t;
         }
      }
      // printf("n0orn1is1: indexPad[%d]=(%d %d %d) t=%f\n",ip,indexPad[ip][irow][0],
      // indexPad[ip][irow][1],indexPad[ip][irow][2],tPads[ip]);
   }
}
