/// \class TpcAlignment
///
/// \brief MiniMC for TPC alignment tests
///
/// \author Valentin Kuzmin, SINP of Moscow State University

#include <iostream>
using namespace std;

#include <TSystem.h>
#include <Riostream.h>
#include "TpcAlignment.h"

void TpcAlignment::Init(BaseTpcSectorGeo &tpcGeo, Int_t Event1, Int_t Event2, Int_t minHits, Int_t *Sectors,
                        const char *InDataFile, const char *aFile)
{
   InHitsFile = InDataFile;
   aInFile    = aFile;
   minNhits   = TMath::Abs(minHits);
   event1     = Event1;
   event2     = Event2;
   fTpcSecGeo = dynamic_cast<TpcSectorGeoAlignmentVAK *>(&tpcGeo); // Null alignment
   if (!fTpcSecGeo) Fatal("TpcAlignment::Init", " !!! Wrong geometry type !!! ");
   if (!aInFile.EqualTo(" ")) { // read a special alignment
      alignmentF    = new TFile(aInFile, "READ");
      tin_alignment = (TTree *)alignmentF->Get("alignment");
      tin_alignment->SetBranchAddress("R0_A", &R0_A);
      tin_alignment->SetBranchAddress("alpha_A", &alpha_A);
      tin_alignment->SetBranchAddress("beta_A", &beta_A);
      tin_alignment->SetBranchAddress("gamma_A", &gamma_A);
      tin_alignment->GetEvent(0);
      //    alignmentF->Close();
      delete tin_alignment;
      delete alignmentF;
      InputHitsFile = new TFile(InHitsFile.Data(), "READ");
   } else {
      InputHitsFile = new TFile(InHitsFile.Data(), "READ");
      tin_alignment = (TTree *)InputHitsFile->Get("alignment");
      tin_alignment->SetBranchAddress("R0_A", &R0_A);
      tin_alignment->SetBranchAddress("alpha_A", &alpha_A);
      tin_alignment->SetBranchAddress("beta_A", &beta_A);
      tin_alignment->SetBranchAddress("gamma_A", &gamma_A);
      tin_alignment->GetEvent(0);
      delete tin_alignment;
   }
   // initialize the TPC Sector geometry with the input alignment
   fTpcSecGeo->TpcSectorGeoA(R0_A, alpha_A, beta_A, gamma_A);

   // read parametrs from the input file
   tin_parameters = (TTree *)InputHitsFile->Get("parameters");
   tin_parameters->SetBranchAddress("b_mf", &b_mf);
   tin_parameters->GetEvent(0);
   delete tin_parameters;

   if ((b_mf[0] * b_mf[0] + b_mf[1] * b_mf[1] + b_mf[2] * b_mf[2]) > 0) {
      tLine = false;
   } else {
      tLine = true;
   }
   // initialize the TPC allignment partial derivetives class
   if (tLine) {
      fLParDer = new DerLMpdTpc(*fTpcSecGeo);
   } else {
      fHParDer = new DerHMpdTpc(*fTpcSecGeo);
   }
   // data about sectors hits
   t_fullChi2 = (TTree *)InputHitsFile->Get("Chi2");
   t_fullChi2->SetBranchAddress("fNhits", &fNhits);
   t_fullChi2->SetBranchAddress("fSumChi2", &fSumChi2);
   t_fullChi2->SetBranchAddress("secNhits", &fNhitSec);
   t_fullChi2->SetBranchAddress("secSumChi2", &fSecChi2);
   t_fullChi2->GetEvent(0);
   delete t_fullChi2;
   // data TRee information from the input file
   tin_data = (TTree *)InputHitsFile->Get("data");

   tin_data->SetBranchAddress("trackPin", &trackPin);
   tin_data->SetBranchAddress("trackPout", &trackPout);
   tin_data->SetBranchAddress("Chi2Fit", &Chi2Fit);
   tin_data->SetBranchAddress("Par0R", &trackDiff[0]);
   tin_data->SetBranchAddress("Par1R", &trackDiff[1]);
   tin_data->SetBranchAddress("Par2R", &trackDiff[2]);
   tin_data->SetBranchAddress("Par3R", &trackDiff[3]);
   tin_data->SetBranchAddress("Par4R", &trackDiff[4]);
   tin_data->SetBranchAddress("Par5R", &trackDiff[5]);
   tin_data->SetBranchAddress("MpdTpcMiniMC_Nhits", &RecoMiniMC_Nhits);
   tin_data->SetBranchAddress("MpdTpcMiniMC_iSec", &RecoMiniMC_iSec);
   tin_data->SetBranchAddress("MpdTpcMiniMC_XHits", &RecoMiniMC_XHits);
   tin_data->SetBranchAddress("MpdTpcMiniMC_YHits", &RecoMiniMC_YHits);
   tin_data->SetBranchAddress("MpdTpcMiniMC_ZHits", &RecoMiniMC_ZHits);
   tin_data->SetBranchAddress("MpdTpcMiniMC_WHits", &RecoMiniMC_WHits);
   tin_data->SetBranchAddress("MpdTpcMiniMC_chi2s", &RecoMiniMC_chi2s);
   //
   for (Int_t i = 0; i < NV; i++) {
      C[i] = 0;
      for (Int_t j = 0; j < NV; j++) A[i][j] = 0;
   }
   for (Int_t i = 0; i < NGV; i++) DFDP[i] = 0;
   for (Int_t i = 0; i < NLV; i++) DFDQ[i] = 0;
   Int_t k = 0;
   for (Int_t is = 0; is < NSEC; is++) {
      for (Int_t i = 0; i < 6; i++) {
         ldp[k] = 10 * is + i;
         k++;
      }
   }
   Int_t spos = 0;
   N          = NLV;
   for (Int_t is = 0; is < 24; is++) {
      sectors[is] = Sectors[is];
      if (sectors[is]) {
         N += 6;
         secpos[is]   = spos;
         possec[spos] = is;
         spos++;
      }
   }
   nS = spos;

   fNhits    = 0; // full hits number
   fSumChi2  = 0; // full chi2 sum of all hits
   proEvents = 0; // proceeded events number
}

//**********************************************************************
//======================================================================
//**********************************************************************

// calculate A & C for the equation
// A*deltaP = C (minimum conditon for the global chi2)
void TpcAlignment::alignment()
{
   if (printL > 0) {
      cout << "miniTpcAlignment: InHitsFile ========  " << InHitsFile.Data() << endl;
      if (!aInFile.EqualTo(" ")) // read a special alignment
         cout << "miniTpcAlignment: TPC alignment from: " << aInFile.Data() << endl;
      else
         cout << "miniTpcAlignment: Using the TPC alignment from the input file: \n" << InHitsFile.Data() << endl;
   }
   int      NLC, NGL = 6, label[6];
   float    derLc[6], derGl[6];
   float    rMeas, sigma;
   Long64_t nhts[NSEC], Nhts = 0;
   for (Int_t i = 0; i < NSEC; i++) nhts[i] = 0;

   // ****************  cicle over events ********************
   // cout<<"-alignment start-"<<endl;
   if (event2 == 0) { // if not the constructor #2
      event1 = 0;
      event2 = tin_data->GetEntries();
   }
   if (printL > 0) cout << "miniTpcAlignment: analyze events: " << event1 << "-" << event2 - 1 << endl;

   if (tLine) {
      NLC = 5;
   } else {
      NLC = 6;
   }

   InputHitsFile->cd();
   for (Int_t iEvt = event1; iEvt < event2; iEvt++) {
      CurrentEvent = iEvt;
      tin_data->GetEvent(iEvt);
      if (tLine) {
         fLParDer->SetLocPar(trackPout);
         TVector3 e  = TVector3(TMath::Sin(trackPout[3]) * TMath::Cos(trackPout[4]),
                                TMath::Sin(trackPout[3]) * TMath::Sin(trackPout[4]), TMath::Cos(trackPout[3]));
         TVector3 x0 = TVector3(trackPout[0], trackPout[1], trackPout[2]);
      } else {
         fHParDer->SetLocPar(trackPout);
      }
      Int_t yes = 0;
      /*
      printf("event=%d nhits=%d  ===============================================\n",
      iEvt,RecoMiniMC_Nhits);
      */
      if (RecoMiniMC_Nhits >= minNhits) {
         for (Int_t ih = 0; ih < RecoMiniMC_Nhits; ih++) {
            Int_t is, spos;
            if (sectors[RecoMiniMC_iSec[ih]]) {
               yes = 1;
               is  = RecoMiniMC_iSec[ih];
            } else
               continue;
            spos = secpos[is];
            nhts[is]++;
            Double_t hit[3] = {RecoMiniMC_XHits[ih], RecoMiniMC_YHits[ih], RecoMiniMC_ZHits[ih]};
            Double_t weight = 1 / (RecoMiniMC_WHits[ih] * RecoMiniMC_WHits[ih]);
            if (tLine) { // line part
                         // cout<<"----------------            // line part"<<endl;
               fLParDer->SetHitTrack(is, hit);
               fLParDer->GetLocDer(dFdq);
               fLParDer->GetGloDer(is, dFdp);
               fLParDer->Getd2Fdpdp(is, d2Fdpdp);
               fLParDer->GetLocDer2(d2Fdqdq);
               fLParDer->GetLocGloDer2(is, d2Fdqdp);
            } else { // helix part
                     // cout<<"----------------            // helix part"<<endl;
               fHParDer->SetHitTrack(is, hit);
               fHParDer->GetLocDer(dFdq);
               fHParDer->GetGloDer(is, dFdp);
               fHParDer->Getd2Fdpdp(is, d2Fdpdp);
               fHParDer->GetLocDer2(d2Fdqdq);
               fHParDer->GetLocGloDer2(is, d2Fdqdp);
            }

            //  fill A & C matrixes

            Int_t l = NLV + 6 * spos;
            for (Int_t i = 0; i < NLV; i++) { // local variables equations
               DFDQ[i] += dFdq[i] * weight;
               C[i] -= dFdq[i] * weight;
               for (Int_t j = 0; j < NLV; j++) {
                  A[i][j] += d2Fdqdq[i][j] * weight;
               }
               for (Int_t j = 0; j < 6; j++) {
                  A[i][l + j] += d2Fdqdp[i][j] * weight;
               }
            }
            for (Int_t i = 0; i < 6; i++) { // global variables equations
               DFDP[l + i] += dFdp[i] * weight;
               C[l + i] -= dFdp[i] * weight;
               for (Int_t j = 0; j < NLV; j++) {
                  A[l + i][j] += d2Fdqdp[j][i] * weight;
               }
               for (Int_t j = 0; j < 6; j++) {
                  A[l + i][l + j] += d2Fdpdp[i][j] * weight;
               }
            }
         }
         if (yes) proEvents++;
         if (printL >= 1) {
            if (proEvents % printL == 0)
               cout << "miniTpcAlignment: events read=" << iEvt << ",   proceeded=" << proEvents << endl;
         }
      }
   }

   //  InputHitsFile->Close();
   delete tin_data;
   delete InputHitsFile;

   N = NLV;
   for (Int_t s = 0; s < 24; s++)
      if (sectors[s]) Nhts += nhts[s];
   for (Int_t s = 0; s < 24; s++)
      if (sectors[s]) N += 6;
   for (Int_t i = 0; i < NLV; i++) B[i] = C[i] / Nhts;

   for (Int_t s = 0; s < 24; s++) {
      if (sectors[s]) {
         Int_t is = secpos[s];
         if (nhts[s]) {
            for (Int_t i = 0; i < 6; i++) {
               B[NLV + 6 * is + i] = C[NLV + 6 * is + i] / Nhts;
               for (Int_t j = 0; j < NLV; j++) A[NLV + 6 * is + i][j] /= Nhts;
               for (Int_t j = NLV + 6 * is; j < NLV + 6 * (s + 1); j++) A[NLV + 6 * is + i][j] /= Nhts;
               for (Int_t j = NLV + 6 * (is + 1); j < NV; j++) A[NLV + 6 * is + i][j] /= Nhts;
            }
         } else
            for (Int_t i = 0; i < 6; i++) {
               for (Int_t j = 0; j < NLV; j++) A[NLV + 6 * is + i][j] = 0;
               for (Int_t j = NLV + 6 * is; j < NLV + 6 * (s + 1); j++)
                  A[NLV + 6 * is + i][j] = (i + 1 + j) * (i + 1 + j);
               for (Int_t j = NLV + 6 * (is + 1); j < NV; j++) A[NLV + 6 * is + i][j] = 0;
               B[NLV + 6 * is + i] = 0;
            }
      }
   }

   if (printL > 1) {
      cout << "==================== A&C ======== N=" << N << " ============" << endl;
      for (Int_t i = 0; i < NLV; i++) {
         for (Int_t m = 0; m < N; m++) printf("%18.10e ", A[i][m]);
         printf("=%18.10e\n", B[i]);
      }
      for (Int_t i = NLV; i < N; i++) {
         Int_t sec;
         if ((i - NLV) % 6 == 0 || i == NLV) {
            sec = (i - NLV) / 6;
            printf("sector=%d   hits=%lld\n", sec, nhts[possec[(i - NLV) / 6]]);
         }
         for (Int_t m = 0; m < NLV; m++) printf("%18.10e ", A[i][m]);
         printf(" ");
         Int_t m1 = NLV + 6 * sec;
         Int_t m2 = m1 + 6;
         m1       = NLV;
         m2       = N;
         for (Int_t m = m1; m < m2; m++) printf("%18.10e ", A[i][m]);
         printf("=%18.10e\n", B[i]);
      }
   }
   linsys(N, X, A, B);
   if (printL > 0) {
      cout << "================= X =================== " << endl;
      if (NLV > 0) printf("q[1-3] (%10.3e,%10.3e,%10.3e) q[4,5] (%10.3e,%10.3e)\n", X[0], X[1], X[2], X[3], X[4]);
      for (Int_t is = 0; is < 24; is++) {
         if (sectors[is]) {
            Int_t ps = secpos[is];
            printf("sector=%2d  ", is);
            Int_t pi = NLV + 6 * ps;
            printf("dR0(%11.5e,%11.5e,%11.5e)  dA(%11.5e,%11.5e,%11.5e)\n", X[pi], X[pi + 1], X[pi + 2], X[pi + 3],
                   X[pi + 4], X[pi + 5]);
         }
      }
      cout << "============== miniTpcAlignment: " << proEvents << " events proceeded ==============" << endl;
   }
}

//**********************************************************************
//======================================================================
//**********************************************************************

//
Int_t TpcAlignment::linsys(Int_t n, Double_t *x, Double_t a[NV][NV], Double_t *b)
{
   //  solves the equation a*x=b

   for (Int_t i = 0; i < n; i++) {
      if (a[i][i] == 0) { // look for string with the non-zero element
         Bool_t d = false;
         for (Int_t j = i + 1; j < n; j++) {
            if (a[j][i] != 0) {                // the i-th diagonal element is non zero
               for (Int_t k = i; k < n; k++) { // change i-th string with j-th row
                  Double_t w = a[j][k];
                  a[j][k]    = a[i][k];
                  a[i][k]    = w;
               }
               Double_t w = b[j]; // change i-th & j-th vector b elements
               b[j]       = b[i];
               b[i]       = w;
               d          = true;
               break; // finish with the i-th string
            }
         }
         if (!d) {
            cout << "linsys_gauss: the degenerate matrix" << endl;
            for (Int_t k = 0; k < n; k++)
               for (Int_t l = 0; l < n; l++)
                  if (l < n - 1)
                     cout << a[k][l] << " ";
                  else
                     cout << a[k][l] << "  " << b[k] << endl;
            return 1;
         }
      }
      // transform to zero the column under the dioganal element
      for (Int_t j = i + 1; j < n; j++) {
         Double_t c = a[j][i] / a[i][i];
         for (Int_t k = i; k < n; k++) {
            a[j][k] -= c * a[i][k];
         }
         b[j] = b[j] - c * b[i];
      }
   }
   // build the solution in x
   for (Int_t i = n - 1; i >= 0; i--) {
      x[i] = b[i];
      for (Int_t j = i + 1; j < n; j++) x[i] -= a[i][j] * x[j];
      x[i] /= a[i][i];
   }
   return 0;
}

//**********************************************************************
//======================================================================
//**********************************************************************

void TpcAlignment::results(Double_t *R0x, Double_t *R0y, Double_t *R0z, Double_t *alpha, Double_t *beta,
                           Double_t *gamma)
{
   for (Int_t s = 0; s < 24; s++)
      if (sectors[s]) {
         R0x[s]   = X[NLV + 6 * secpos[s]];
         R0y[s]   = X[NLV + 6 * secpos[s] + 1];
         R0z[s]   = X[NLV + 6 * secpos[s] + 2];
         alpha[s] = X[NLV + 6 * secpos[s] + 3];
         beta[s]  = X[NLV + 6 * secpos[s] + 4];
         gamma[s] = X[NLV + 6 * secpos[s] + 5];
      } else {
         R0x[s]   = 0;
         R0y[s]   = 0;
         R0z[s]   = 0;
         alpha[s] = 0;
         beta[s]  = 0;
         gamma[s] = 0;
      }
}

//**********************************************************************
//======================================================================
//**********************************************************************

void TpcAlignment::getDFDP(Double_t *der)
{
   for (Int_t i = 0; i < NLV; i++) der[i] = DFDQ[i];
   // for(Int_t s=0;s<nS;s++)
   // cout<<s<<" "<<possec[i]<<"   "<<NLV+6*possec[i]<<" "<<NLV+6*s+i<<endl;
   //   for(Int_t i=0;i<6;i++)
   //    der[NLV+6*possec[s]+i]=DFDP[NLV+6*s+i];
   for (Int_t s = 0; s < 24; s++)
      if (sectors[s]) {
         for (Int_t i = 0; i < 6; i++) der[NLV + 6 * s + i] = DFDP[NLV + 6 * secpos[s] + i];
      }
}

//**********************************************************************
//======================================================================
//**********************************************************************

// scan events, pass them to Millepede calculate A & C for the equation
// A*deltaP = C (minimum conditon for the global chi2)
void TpcAlignment::Alignment(Double_t &F, Double_t *gradF)
{
   if (printL > 0) {
      cout << "miniTpcAlignment: InHitsFile ========  " << InHitsFile.Data() << endl;
      if (!aInFile.EqualTo(" ")) // read a special alignment
         cout << "miniTpcAlignment: TPC alignment from: " << aInFile.Data() << endl;
      else
         cout << "miniTpcAlignment: Using the TPC alignment from the input file: \n" << InHitsFile.Data() << endl;
   }
   int      NLC, NGL       = 6, label[6];
   Long64_t nhts[24], Nhts = 0;
   for (Int_t i = 0; i < 24; i++) nhts[i] = 0;
   for (Int_t i = 0; i < NLV; i++) DFDQ[i] = 0;
   for (Int_t i = 0; i < NGV; i++) DFDP[i] = 0;
   Double_t Fun = 0;

   // ****************  cicle over events ********************
   if (event2 == 0) { // if not the constructor #2
      event1 = 0;
      event2 = tin_data->GetEntries();
   }
   if (printL > 0) cout << "miniTpcAlignment: analyze events: " << event1 << "-" << event2 - 1 << endl;
   if (tLine) {
      NLC = 5;
   } else {
      NLC = 6;
   }
   InputHitsFile->cd();
   for (Int_t iEvt = 0; iEvt < event1; iEvt++) tin_data->GetEvent(iEvt);
   for (Int_t iEvt = event1; iEvt < event2; iEvt++) {
      tin_data->GetEvent(iEvt);
      fLParDer->SetLocPar(trackPout);
      Int_t    yes = 0;
      TVector3 x1, res;
      TVector3 e  = TVector3(TMath::Sin(trackPout[3]) * TMath::Cos(trackPout[4]),
                             TMath::Sin(trackPout[3]) * TMath::Sin(trackPout[4]), TMath::Cos(trackPout[3]));
      TVector3 x0 = TVector3(trackPout[0], trackPout[1], trackPout[2]);
      // printf("---e=%d nh=%d\n",iEvt,RecoMiniMC_Nhits);
      if (RecoMiniMC_Nhits >= minNhits) {
         for (Int_t ih = 0; ih < RecoMiniMC_Nhits; ih++) {
            Int_t is, spos;
            if (sectors[RecoMiniMC_iSec[ih]]) {
               yes = 1;
               is  = RecoMiniMC_iSec[ih];
            } else
               continue;
            spos = secpos[is];
            nhts[is]++;
            Double_t hit[3] = {RecoMiniMC_XHits[ih], RecoMiniMC_YHits[ih], RecoMiniMC_ZHits[ih]};
            x1              = TVector3(hit[0], hit[1], hit[2]);
            res             = e.Cross(x0 - x1);
            Double_t sigma  = RecoMiniMC_WHits[ih];
            Double_t delta  = res.Mag() / sigma;
            Double_t weight = 1 / (sigma * sigma);
            //        weight=1;
            Fun += delta * delta;
            if (tLine) {
               fLParDer->SetHitTrack(is, hit);
               fLParDer->GetLocDer(dFdq);
               //          fLParDer->GetGloDerT(is,dFdp);
               fLParDer->GetGloDer(is, dFdp);
               // if(10==10&&iEvt==event1)
               // if(is==1)printf("s=%2d dFdp(%g,%g,%g,%g,%g,%g)\n",is,dFdp[0],dFdp[1],dFdp[2],dFdp[3],dFdp[4],dFdp[5]);
            } else {
               // helix part
            }

            Int_t l = NLV + 6 * spos;
            for (Int_t i = 0; i < NLV; i++) { // local variables equations
               DFDQ[i] += dFdq[i] * weight;
            }
            for (Int_t i = 0; i < 6; i++) { // global variables equations
               DFDP[l + i] += dFdp[i] * weight;
            }
            // if(10==10&&iEvt==event1)
            // if(is==1)printf("DFDP(%g,%g,%g,%g,%g,%g)\n\n",/DFDP[l+0],DFDP[l+1],DFDP[l+2],DFDP[l+3],DFDP[l+4],DFDP[l+5]);
         }
         if (yes) proEvents++;
         if (printL > 1) {
            if (proEvents % printL == 0)
               cout << "miniTpcAlignment: events read=" << iEvt << ",   proceeded=" << proEvents << endl;
         }
      }
      // if(iEvt==event1+10)for(Int_t i=0;i<NGV;i++) printf("DFDP[%d]=%g\n",i,DFDP[i]);
   }
   // for(Int_t i=0;i<NGV;i++) printf("DFDP[%d]=%g\n",i,DFDP[i]);

   delete tin_data;
   delete InputHitsFile;

   N = NLV;
   for (Int_t s = 0; s < 24; s++) {
      if (sectors[s]) {
         Nhts += nhts[s];
         N += 6;
      }
   }
   Double_t glength = 0;
   for (Int_t i = 0; i < NLV; i++) glength += DFDQ[i] * DFDQ[i];
   for (Int_t i = 0; i < NGV; i++) glength += DFDP[i] * DFDP[i];
   //  glength=TMath::Sqrt(glength);
   glength = Nhts;
   for (Int_t i = 0; i < NLV; i++) gradF[i] = DFDQ[i] / glength;
   for (Int_t i = 0; i < NGV; i++) gradF[i + NLV] = DFDP[i] / glength;
   F = Fun / Nhts;
   if (printL > 1) {
      cout << "miniTpcAlignment: ====== dFdq,dFdp ====== N=" << N << "(" << NLV << "+" << NGV << ") =====" << endl;
      for (Int_t i = 0; i < NLV; i++) {
         printf("dFdq[%d]=%18.10e\n", i + 1, DFDQ[i]);
      }
      for (Int_t i = 0; i < NGV; i++) {
         Int_t sec;
         if (i % 6 == 0) {
            sec = i / 6;
            printf("sector=%d   hits=%lld\n", sec, nhts[possec[sec]]);
            for (Int_t m = 0; m < 6; m++) {
               if (m < 3) printf("x[%d]=%18.10e  ", m, R0_A[m][sec]);
               if (m == 3) printf("x[%d]=%18.10e  ", m, alpha_A[sec]);
               if (m == 4) printf("x[%d]=%18.10e  ", m, beta_A[sec]);
               if (m == 5) printf("x[%d]=%18.10e  ", m, gamma_A[sec]);
               printf("DFDP[%d]=%18.10e\n", m, DFDP[i + m]);
            }
         }
      }
      printf("miniTpcAlignment: F/Nhits=%18.15f\n", F);
   }
   //  for(Int_t i=0;i<N;i++) B[i]=C[i];
   if (printL > 0) {
      cout << "============== miniTpcAlignment: " << proEvents << " events proceeded ==============" << endl;
   }
}
