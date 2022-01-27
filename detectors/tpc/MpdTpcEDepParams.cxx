/// \class MpdTpcEDepParams
///
/// Parameters for MPD TPC dE/dx simulation
/// \author Igor Rufanov (LHEP, JINR, Dubna) - development
/// \author Alexander Zinchenko (LHEP, JINR, Dubna) - porting to MpdRoot
/// 18/11/2019

#include "MpdTpcEDepParams.h"

#include <TFile.h>
#include <TH2D.h>
#include <TMath.h>
#include <TRandom.h>
//#include <TSystem.h>

#include <iostream>

using std::cout;
using std::endl;
using std::vector;

MpdTpcEDepParams *MpdTpcEDepParams::fgTpcEdepParams = nullptr;

//__________________________________________________________________________

MpdTpcEDepParams *MpdTpcEDepParams::Instance()
{
   /// Get pointer to the TPC ionization loss patameter object

   if (!fgTpcEdepParams) {
      fgTpcEdepParams = new MpdTpcEDepParams;

      // Automatic deletion
      std::atexit(DestroyInstance);
   }
   return fgTpcEdepParams;
}

//__________________________________________________________________________

void MpdTpcEDepParams::Init()
{
   /// Get histograms and initialize parameters

   TFile *file = new TFile("$VMCWORKDIR/input/TpcEDepParamsHeed.root");
   TH1D * h    = (TH1D *)file->Get("NCollTab");

   fNCollDensBgBins  = h->GetNbinsX();
   fCollDensL10BgMin = h->GetBinCenter(1);
   // fCollDensL10BgMin = -0.5;
   fCollDensL10BgMax  = h->GetBinCenter(fNCollDensBgBins);
   fCollDensL10BgStep = (fCollDensL10BgMax - fCollDensL10BgMin) / (fNCollDensBgBins - 1);
   fCollDens          = new Double_t[fNCollDensBgBins];

   for (Int_t i = 0; i < fNCollDensBgBins; i++) fCollDens[i] = h->GetBinContent(i + 1);

   TH2D *he      = (TH2D *)file->Get("InvertedProbFuncEDep2");
   fNEneBgBins   = he->GetNbinsX();
   fEneL10BgMin  = he->GetXaxis()->GetBinCenter(1);
   fEneL10BgMax  = he->GetXaxis()->GetBinCenter(fNEneBgBins);
   fEneL10BgStep = (fEneL10BgMax - fEneL10BgMin) / (fNEneBgBins - 1);

   fNEneProbBins = he->GetNbinsY();
   fEneProbStep  = 1.0 / (fNEneProbBins - 2);

   vector<Double_t> probs(fNEneProbBins);
   fEne.resize(fNEneBgBins, probs);

   for (Int_t i = 0; i < fNEneBgBins; i++) {
      for (Int_t j = 0; j < fNEneProbBins; j++) {
         // fEne[i][j] = he->GetBinContent(i+1,j+1);
         fEne[i][j] = he->GetBinContent(i + 1, j + 2);
      }
   }
   file->Close();
}

//__________________________________________________________________________

Double_t MpdTpcEDepParams::GetCollisionDensity(Double_t log10bg)
{
   /// Number of clusters / cm vs log10(beta*gamma)

   Double_t bin = (log10bg - fCollDensL10BgMin) / fCollDensL10BgStep;
   if (bin <= 0.)
      return fCollDens[0];
   else if (bin >= fNCollDensBgBins - 2)
      return fCollDens[fNCollDensBgBins - 2];
   else {
      Int_t    ib = (Int_t)bin;
      Double_t d  = bin - ib;
      return fCollDens[ib] * (1.0 - d) + fCollDens[ib + 1] * d;
   }
}

//__________________________________________________________________________

Double_t MpdTpcEDepParams::GetRandEnergy(Double_t log10bg)
{
   /// Energy in the cluster vs log10(beta*gamma)

   Int_t    ib  = 0; //, ib1 = 0;
   Double_t bin = (log10bg - fEneL10BgMin) / fEneL10BgStep;
   if (bin <= 0.)
      ib = 0;
   else if (bin >= fNEneBgBins - 1)
      ib = fNEneBgBins - 1;
   else {
      // ib = round( bin);
      // smearing at bin boundary
      ib         = (Int_t)bin;
      Double_t d = bin - ib;
      if (gRandom->Uniform(0., 1.) > (1.0 - d)) ib++;
      if (ib >= fNEneBgBins) exit(0);
   }
   /*
   else {
     ib = (Int_t) (bin + 0.5);
     if (ib >= fNEneBgBins) exit(0);
   }
   */
   /*
   else {
     // Linear combination
     ib = ib1 = (Int_t) bin;
     if (ib1 < fNEneBgBins - 1) ++ib1;
   }
   */

   Double_t prob    = gRandom->Uniform(0., 1.);
   Double_t probBin = prob / fEneProbStep;
   Int_t    ip      = (Int_t)probBin;
   Double_t p       = probBin - ip;

   Double_t l10 = fEne[ib][ip] * (1. - p) + fEne[ib][ip + 1] * p;
   // Double_t l101 = fEne[ib1][ip] * (1.-p) + fEne[ib1][ip+1] * p;
   // Double_t d = bin - ib;
   // l10 = l10 * (1.0 - d) + l101 * d;
   if (ip + 1 >= fNEneProbBins) exit(0);
   // cout << l10 << " " << ib << " " << ip << endl;
   return TMath::Power(10., l10);

   /*
   Double_t l10 = TMath::Power(10., fEne[ib][ip]);
   Double_t l101 = TMath::Power(10., fEne[ib][ip+1]);
   return l10 * (1.0 - p) + l101 * p;
   */
}

//__________________________________________________________________________

Double_t MpdTpcEDepParams::GetEloss(Double_t log10bg, Double_t charge, Double_t step)
{
   /// Energy loss for given log10(beta*gamma) and step

   Double_t collDens = GetCollisionDensity(log10bg) * charge * charge;
   Double_t nColl    = gRandom->PoissonD(collDens * step);
   Double_t eloss    = 0;

   for (Int_t i = 0; i < nColl; ++i) eloss += GetRandEnergy(log10bg);
   // cout << nColl << " xxxx " << eloss << endl;
   return eloss * 1e-9;
}

//__________________________________________________________________________