//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class ParticleIdentification
//      see ParticleIdentification.h for details
//
// Environment:
//      Software developed for MPD at NICA.
//
// Author List:
//      Gyulnara Eyyubova
//
//-----------------------------------------------------------

// C/C++ Headers ----------------------
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>
// This Class' Header ------------------
#include "MpdTOFpid.h"

using namespace std;
using namespace TMath;

// Class Member definitions -----------

MpdTOFpid::MpdTOFpid()
   : fnBeta(1.), ProtonParHyper(), PionParHyper(), KaonParHyper(), ElectronParHyper(), sigmasTofPi(), sigmasTofPr(),
     sigmasTofKa(), sigmasTofEl()
{
   ReadTOFResponse();
}

MpdTOFpid::~MpdTOFpid() {}

Int_t MpdTOFpid::GetTofProbs(Float_t P, Float_t beta, Float_t &Ppi, Float_t &Pk, Float_t &Pp, Float_t &Pe, Int_t method)
{

   const Int_t Ntypes = 4; // order: pion, kaon, proton, electron
   Float_t     resultProb[Ntypes];

   const Int_t Nintervals            = 72;
   Double_t    P_int[Nintervals + 1] = {
      0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.3,  0.33, 0.36, 0.39, 0.42, 0.45, 0.48, 0.51, 0.54, 0.57, 0.6,  0.63,
      0.66, 0.69, 0.72, 0.75, 0.78, 0.81, 0.84, 0.87, 0.9,  0.93, 0.96, 0.99, 1.02, 1.05, 1.08, 1.11, 1.14, 1.17, 1.2,
      1.23, 1.26, 1.29, 1.32, 1.35, 1.38, 1.41, 1.44, 1.47, 1.5,  1.53, 1.56, 1.59, 1.62, 1.65, 1.68, 1.71, 1.74, 1.77,
      1.8,  1.83, 1.86, 1.89, 1.92, 1.95, 1.98, 2.01, 2.04, 2.07, 2.1,  2.13, 2.16, 2.19, 2.22, 3};

   if (Abs(beta) < 1e-6) return 1;
   const Float_t oneOverBeta = 1.0 / beta;

   Float_t sum = 0.0; // for normalizing

   /*A priori coefficients for Bayes approach */
   Float_t bayesAprioriCoefficients[Ntypes];

   /*A "measured" probabilities */
   Float_t gausProb[Ntypes];
   for (Int_t i = 0; i < Ntypes; i++) gausProb[i] = 0;

   Int_t i_p = 0; // momentum interval

   switch (method) {

   case 0: // Hyperbolic (1/beta vs p) with bayesian coefficients

      bayesAprioriCoefficients[0] = 1.0;
      bayesAprioriCoefficients[1] = 1.0;
      bayesAprioriCoefficients[2] = 1.0;
      bayesAprioriCoefficients[3] = 1.0;

      if (P > 3) i_p = 71;
      for (Int_t k = 0; k < Nintervals; k++)
         if (P >= P_int[k] && P < P_int[k + 1]) i_p = k;
      gausProb[0] = Gaus(oneOverBeta, HyperbolicFunction(P, PionParHyper), sigmasTofPi[i_p], kTRUE);
      gausProb[1] = Gaus(oneOverBeta, HyperbolicFunction(P, KaonParHyper), sigmasTofKa[i_p], kTRUE);
      gausProb[2] = Gaus(oneOverBeta, HyperbolicFunction(P, ProtonParHyper), sigmasTofPr[i_p], kTRUE);
      gausProb[3] = Gaus(oneOverBeta, HyperbolicFunction(P, ElectronParHyper), sigmasTofEl[i_p], kTRUE);

      sum = gausProb[0] + gausProb[1] + gausProb[2] + gausProb[3];
      if (sum == 0.) return 1;
      gausProb[0] /= sum;
      gausProb[1] /= sum;
      gausProb[2] /= sum;
      gausProb[3] /= sum;
      break;

   case 1: // Hyperbolic (1/beta vs p) with special byesian coefficients

      if (P > 3) i_p = 71;
      for (Int_t k = 0; k < Nintervals; k++)
         if (P >= P_int[k] && P < P_int[k + 1]) i_p = k;

      if (P < 0.3) {
         bayesAprioriCoefficients[0] = 0.28;
         bayesAprioriCoefficients[1] = 0.25;
         bayesAprioriCoefficients[2] = 0.26;
         bayesAprioriCoefficients[3] = 0.21;
      } else if ((P >= 0.3) && (P < 0.4)) {
         bayesAprioriCoefficients[0] = 0.91;
         bayesAprioriCoefficients[1] = 0.02;
         bayesAprioriCoefficients[2] = 0.06;
         bayesAprioriCoefficients[3] = 0.01;
      } else if ((P >= 0.4) && (P < 0.5)) {
         bayesAprioriCoefficients[0] = 0.70;
         bayesAprioriCoefficients[1] = 0.11;
         bayesAprioriCoefficients[2] = 0.13;
         bayesAprioriCoefficients[3] = 0.06;
      } else if ((P >= 0.5) && (P < 0.6)) {
         bayesAprioriCoefficients[0] = 0.39;
         bayesAprioriCoefficients[1] = 0.19;
         bayesAprioriCoefficients[2] = 0.32;
         bayesAprioriCoefficients[3] = 0.10;
      } else if ((P >= 0.6) && (P < 0.7)) {
         bayesAprioriCoefficients[0] = 0.41;
         bayesAprioriCoefficients[1] = 0.17;
         bayesAprioriCoefficients[2] = 0.37;
         bayesAprioriCoefficients[3] = 0.5;
      } else if ((P >= 0.7) && (P < 0.8)) {
         bayesAprioriCoefficients[0] = 0.43;
         bayesAprioriCoefficients[1] = 0.16;
         bayesAprioriCoefficients[2] = 0.31;
         bayesAprioriCoefficients[3] = 0.10;
      } else if ((P >= 0.8) && (P < 0.9)) {
         bayesAprioriCoefficients[0] = 0.44;
         bayesAprioriCoefficients[1] = 0.11;
         bayesAprioriCoefficients[2] = 0.43;
         bayesAprioriCoefficients[3] = 0.02;
      } else if ((P >= 0.9) && (P < 1.0)) {
         bayesAprioriCoefficients[0] = 0.36;
         bayesAprioriCoefficients[1] = 0.21;
         bayesAprioriCoefficients[2] = 0.36;
         bayesAprioriCoefficients[3] = 0.07;
      } else if ((P >= 1.0) && (P < 1.2)) {
         bayesAprioriCoefficients[0] = 0.32;
         bayesAprioriCoefficients[1] = 0.32;
         bayesAprioriCoefficients[2] = 0.32;
         bayesAprioriCoefficients[3] = 0.04;
      } else if (P >= 1.2) {
         bayesAprioriCoefficients[0] = 0.34;
         bayesAprioriCoefficients[1] = 0.27;
         bayesAprioriCoefficients[2] = 0.31;
         bayesAprioriCoefficients[3] = 0.08;
      }

      gausProb[0] = Gaus(oneOverBeta, HyperbolicFunction(P, PionParHyper), sigmasTofPi[i_p], kTRUE);
      gausProb[1] = Gaus(oneOverBeta, HyperbolicFunction(P, KaonParHyper), sigmasTofKa[i_p], kTRUE);
      gausProb[2] = Gaus(oneOverBeta, HyperbolicFunction(P, ProtonParHyper), sigmasTofPr[i_p], kTRUE);
      gausProb[3] = Gaus(oneOverBeta, HyperbolicFunction(P, ElectronParHyper), sigmasTofEl[i_p], kTRUE);

      sum = gausProb[0] + gausProb[1] + gausProb[2] + gausProb[3];
      if (sum == 0.) return 1;
      gausProb[0] /= sum;
      gausProb[1] /= sum;
      gausProb[2] /= sum;
      gausProb[3] /= sum;
      break;

   case 2: // n-sigma method,
      Float_t delta[Ntypes];

      bayesAprioriCoefficients[0] = 1.0;
      bayesAprioriCoefficients[1] = 1.0;
      bayesAprioriCoefficients[2] = 1.0;
      bayesAprioriCoefficients[3] = 1.0;

      if (P > 3) i_p = 71;
      for (Int_t k = 0; k < Nintervals; k++)
         if (P >= P_int[k] && P < P_int[k + 1]) i_p = k;

      delta[0] = TMath::Abs(oneOverBeta - HyperbolicFunction(P, PionParHyper)) / sigmasTofPi[i_p];
      delta[1] = TMath::Abs(oneOverBeta - HyperbolicFunction(P, KaonParHyper)) / sigmasTofKa[i_p];
      delta[2] = TMath::Abs(oneOverBeta - HyperbolicFunction(P, ProtonParHyper)) / sigmasTofPr[i_p];
      delta[3] = TMath::Abs(oneOverBeta - HyperbolicFunction(P, ElectronParHyper)) / sigmasTofEl[i_p];

      for (Int_t ii = 0; ii < 4; ii++) {
         if (delta[ii] < fnBeta) gausProb[ii] = 1;
      }
      Int_t c = 0;
      for (Int_t ii = 0; ii < 4; ii++) {
         c = gausProb[ii];
         for (Int_t jj = ii + 1; jj < 4; jj++) {
            if (c == 1 && gausProb[jj] == 1) {
               gausProb[ii] = 0;
               gausProb[jj] = 0;
            }
         }
      }

      break;
   }

   if (method != 0 && method != 1 && method != 2) {
      printf("No TOF method matches.\n");
   }

   else {
      Ppi = 0.;
      Pk  = 0.;
      Pp  = 0.;
      Pe  = 0.;
      if (BayesFunction(gausProb, bayesAprioriCoefficients, resultProb, Ntypes) != 1) {
         Ppi = resultProb[0];
         Pk  = resultProb[1];
         Pp  = resultProb[2];
         Pe  = resultProb[3];
      }
   }
   return 0;
}

Float_t MpdTOFpid::HyperbolicFunction(Float_t x, Float_t *p)
{
   return Sqrt(1 + p[0] * p[0] / (x * x)) + p[1]; // p[0]=m, p[1]=additive const
}

/*Bayes function for probabilities calculation*/
Int_t MpdTOFpid::BayesFunction(Float_t *measProb, Float_t *aprioriProb, Float_t *bayesProb, Int_t N)
{

   /*measProb    - measured probabilities from TOF/TOF/etc*/
   /*aprioriProb - a priori probabilities from some magic*/
   /*bayesProb   - result bayes probabilities */
   /*N           - number of types {pion, kaon, proton, electron} */

   Float_t sumBayes = 0.0;
   for (Int_t i = 0; i < N; ++i) {
      sumBayes += measProb[i] * aprioriProb[i];
   }

   if (sumBayes == 0.) return 1;

   for (Int_t i = 0; i < N; ++i) {
      bayesProb[i] = measProb[i] * aprioriProb[i] / sumBayes;
   }

   return 0;
}

void MpdTOFpid::ReadTOFResponse()
{

   std::string s(100, ' ');

   /* Accessing the oneOverBeta response */
   TString  dir      = getenv("VMCWORKDIR");
   TString  BetaFile = dir + "/input/MPDTOFPidResponseVhelle.txt";
   ifstream input;
   input.open(BetaFile);

   if (!input.is_open()) {
      cerr << "PID: READ TOF RESPONSE: Cannot open file " << BetaFile.Data() << endl;
   }

   // const Int_t Ntypes = 4; //order: pion, kaon, proton, electron

   if (input) {

      getline(input, s);
      for (Int_t j = 0; j < 2; j++) input >> PionParHyper[j];
      getline(input, s);
      getline(input, s);
      for (Int_t j = 0; j < 2; j++) input >> KaonParHyper[j];
      getline(input, s);
      getline(input, s);
      for (Int_t j = 0; j < 2; j++) input >> ProtonParHyper[j];
      getline(input, s);
      getline(input, s);
      for (Int_t j = 0; j < 2; j++) input >> ElectronParHyper[j];

      getline(input, s);
      getline(input, s);
      for (Int_t j = 0; j < 72; j++) input >> sigmasTofPi[j];
      getline(input, s);
      getline(input, s);
      for (Int_t j = 0; j < 72; j++) input >> sigmasTofKa[j];
      getline(input, s);
      getline(input, s);
      for (Int_t j = 0; j < 72; j++) input >> sigmasTofPr[j];
      getline(input, s);
      getline(input, s);
      for (Int_t j = 0; j < 72; j++) input >> sigmasTofEl[j];

      input.close();
   }

   //  Draw response
   /*
       const Int_t Nintervals = 72;
       Double_t P_int[Nintervals + 1] =
   {0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.3,0.33,0.36,0.39,0.42,0.45,0.48,0.51,0.54,0.57,0.6,0.63,0.66,0.69,0.72,0.75,0.78,0.81,0.84,0.87,0.9,0.93,0.96,0.99,1.02,1.05,1.08,
   1.11,1.14,1.17,1.2,1.23,1.26,1.29,1.32,1.35,1.38,1.41,1.44,1.47,1.5,1.53,1.56,1.59,1.62,1.65,1.68,1.71,1.74,1.77,1.8,1.83,1.86,1.89,1.92,1.95,1.98,2.01,2.04,2.07,2.1,
   2.13,2.16,2.19,2.22, 3};

     Float_t P[72], P_err[72], mean_pi[72], mean_pr[72], mean_ka[72], mean_el[72];


     for(Int_t i=0; i<72; i++){
      P[i]=P_int[i]+0.5*(P_int[i+1]-P_int[i]);
      P_err[i]=0;
      mean_pi[i]=HyperbolicFunction(P[i], PionParHyper);
      mean_ka[i]=HyperbolicFunction(P[i], KaonParHyper);
      mean_pr[i]=HyperbolicFunction(P[i], ProtonParHyper);
      mean_el[i]=HyperbolicFunction(P[i], ElectronParHyper);
      }


      TGraphErrors *pi= new TGraphErrors(72, P, mean_pi, P_err, sigmasTofPi);
      pi->SetMarkerColor(2);
      pi->SetLineColor(2);
      pi->SetMarkerStyle(20);
      pi->SetMarkerSize(0.5);
      pi->SetFillColor(1);
      pi->SetFillStyle(3002);

      TGraphErrors *ka= new TGraphErrors(72, P, mean_ka, P_err, sigmasTofKa);
      ka->SetMarkerColor(2);
      ka->SetLineColor(2);
      ka->SetMarkerStyle(20);
      ka->SetMarkerSize(0.5);
      ka->SetFillColor(1);
      ka->SetFillStyle(3002);

      TGraphErrors *pr= new TGraphErrors(72, P, mean_pr, P_err, sigmasTofPr);
      pr->SetMarkerColor(2);
      pr->SetLineColor(2);
      pr->SetMarkerStyle(20);
      pr->SetMarkerSize(0.5);
      pr->SetFillColor(1);
      pr->SetFillStyle(3002);

      TGraphErrors *el= new TGraphErrors(72, P, mean_el, P_err, sigmasTofEl);
      el->SetMarkerColor(2);
      el->SetLineColor(2);
      el->SetMarkerStyle(20);
      el->SetMarkerSize(0.5);
      el->SetFillColor(1);
      el->SetFillStyle(3002);

      pr->Draw("AC3");
      ka->Draw("3Csame");
      pi->Draw("3Csame");
      el->Draw("3Csame");

     // pr->Draw("LXsame");
     // ka->Draw("LXsame");
     // pi->Draw("LXsame");
     // el->Draw("LXsame");

   */
}

ClassImp(MpdTOFpid);
