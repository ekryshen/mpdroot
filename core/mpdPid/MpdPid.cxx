#include "MpdPid.h"

MpdPid::MpdPid() : TObject()
{
   fGaus      = nullptr;
   fAsymGaus  = nullptr;
   fGaus2     = nullptr;
   fAsymGaus2 = nullptr;

   for (Int_t iType = 0; iType < MpdPidUtils::kNSpecies; iType++) {
      fProb[iType]                = 0.0;
      fEnLossSigmasArray[iType]   = 0.0;
      fMSquaredSigmasArray[iType] = 0.0;
      fBayesCoefficients[iType]   = 0.0;
      fNSigSpecies[iType]         = kFALSE;
      fdEdxBottomBound[iType]     = 0.0;
      fdEdxTopBound[iType]        = 0.0;
      fm2BottomBound[iType]       = 0.0;
      fm2TopBound[iType]          = 0.0;
   }
   fMethod        = kTRUE;
   fPrRatio       = 1.0;
   fEnergy        = 4.0;
   fTrackingState = MpdPidUtils::kCFHM;
   fCharge        = MpdPidUtils::kPos;
}

MpdPid::MpdPid(Double_t sigmaTof, Double_t sigmaEloss, Double_t sqrts, Double_t EnLossCoef, TString Generator,
               TString Tracking, TString NSigPart)
   : TObject(), fCharge(MpdPidUtils::kPos), fEnergy(sqrts)
{
   for (Int_t iType = 0; iType < MpdPidUtils::kNSpecies; iType++) {
      fdEdxBottomBound[iType] = sigmaEloss;
      fdEdxTopBound[iType]    = sigmaEloss;
      fm2BottomBound[iType]   = sigmaTof;
      fm2TopBound[iType]      = sigmaTof;
   }
   Init(Generator, Tracking, NSigPart, EnLossCoef);
}

MpdPid::~MpdPid()
{
   for (auto it = fdEdxBBMap.begin(); it != fdEdxBBMap.end(); ++it) {
      for (TF1 *F : it->second) delete F;
   }
   for (auto it = fdEdxSigmaMap.begin(); it != fdEdxSigmaMap.end(); ++it) {
      for (TF1 *F : it->second) delete F;
   }
   for (auto it = fdEdxDeltaMap.begin(); it != fdEdxDeltaMap.end(); ++it) {
      for (TF1 *F : it->second) delete F;
   }
   for (auto it = fParM2Map.begin(); it != fParM2Map.end(); ++it) {
      for (TF1 *F : it->second) delete F;
   }
   for (auto it = fPartYieldMap.begin(); it != fPartYieldMap.end(); ++it) {
      for (TF1 *F : it->second) delete F;
   }
   delete fGaus;
   delete fAsymGaus;
   delete fGaus2;
   delete fAsymGaus2;
}

Double_t MpdPid::parElBB(Double_t *x, Double_t *par)
{
   Double_t ans = par[0] / x[0] + par[1];
   if (fTrackingState == MpdPidUtils::kCF) {
      if (x[0] < 0.15) ans *= 0.980;
      if (x[0] < 0.10) ans *= 1.025;
   }
   return ans;
}

Double_t MpdPid::parMuBB(Double_t *x, Double_t *par)
{
   Double_t x1  = par[0] / TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.011), par[3]);
   Double_t x2  = par[1] - TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.011), par[3]);
   Double_t x3  = TMath::Log(par[2] + TMath::Power(1.0 / (x[0] / 0.1057), par[4]));
   Double_t ans = x1 * (x2 - x3);
   if (fTrackingState == MpdPidUtils::kCF) {
      if (x[0] < 0.10) ans *= 0.960;
      if (x[0] > 1.25) ans *= 1.005;
      if (x[0] > 1.90) ans *= 0.985;
   } else if (fTrackingState == MpdPidUtils::kHP) {
      if (x[0] < 0.15) ans *= 0.970;
      if (x[0] < 0.20) ans *= 1.015;
      if (x[0] > 1.35) ans *= 1.005;
      if (x[0] > 1.60) ans *= 1.005;
      if (x[0] > 1.90) ans *= 0.985;
   } else
      ans *= 1.0;
   return ans;
}

Double_t MpdPid::parPiBB(Double_t *x, Double_t *par)
{
   Double_t x1  = par[0] / TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.01949), par[3]);
   Double_t x2  = par[1] - TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.01949), par[3]);
   Double_t x3  = TMath::Log(par[2] + TMath::Power(1.0 / (x[0] / 0.1396), par[4]));
   Double_t ans = x1 * (x2 - x3);
   if (fTrackingState == MpdPidUtils::kCFHM) {
      // if ( ( x[0] >= 0.05 ) && ( x[0] < 0.10 ) ) ans *= 1.025;
      // if ( ( x[0] >= 0.10 ) && ( x[0] < 0.15 ) ) ans *= 0.970;
      if ((x[0] >= 0.05) && (x[0] < 0.10)) ans *= 0.98;
   } else if (fTrackingState == MpdPidUtils::kCF) {
      if ((x[0] > 0.05) && (x[0] < 0.10)) ans *= 0.830;
      if ((x[0] > 0.10) && (x[0] < 0.15)) ans *= 0.978;
      if ((x[0] > 0.15) && (x[0] < 0.20)) ans *= 1.012;
      if ((x[0] > 0.25) && (x[0] < 0.30)) ans *= 0.987;
      if ((x[0] > 2.60) && (x[0] < 3.00)) ans *= 0.980;
   } else {
      if (x[0] < 0.10) ans *= 0.895;
      if (x[0] < 0.15) ans *= 0.955;
      if (x[0] < 0.25) ans *= 1.020;
      if (x[0] < 0.30) ans *= 0.990;
      if (x[0] > 1.70) ans *= 1.005;
      if (x[0] > 2.50) ans *= 0.990;
   }
   return ans;
}

Double_t MpdPid::parKaBB(Double_t *x, Double_t *par)
{
   Double_t x1, x2, x3, p[5], xint = 0.50933974, ans;
   for (Int_t k = 0; k < 5; k++) p[k] = x[0] < xint ? par[k] : par[k + 5];
   x1  = p[0] / TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.2437), p[3]);
   x2  = p[1] - TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.2437), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / (x[0] / 0.4937), p[4]));
   ans = x1 * (x2 - x3);
   if (fTrackingState == MpdPidUtils::kCFHM) {
      if ((x[0] >= 0.05) && (x[0] < 0.15)) ans *= 0.88;
      if ((x[0] >= 0.15) && (x[0] < 0.20)) ans *= 0.95;
      if ((x[0] >= 0.20) && (x[0] < 0.25)) ans *= 0.99;
      if ((x[0] >= 2.60) && (x[0] < 3.00)) ans *= 0.99;
      // if ( ( x[0] >= 0.05 ) && ( x[0] < 0.10 ) ) ans *= 0.945;
      // if ( ( x[0] >= 0.10 ) && ( x[0] < 0.15 ) ) ans *= 0.930;
      // if ( ( x[0] >= 0.15 ) && ( x[0] < 0.20 ) ) ans *= 0.970;
      // if ( ( x[0] >= 0.20 ) && ( x[0] < 0.25 ) ) ans *= 0.990;
      // if ( ( x[0] >= 1.90 ) && ( x[0] < 2.10 ) ) ans *= 0.990;
      // if ( ( x[0] >= 2.10 ) && ( x[0] < 2.20 ) ) ans *= 0.985;
      // if ( ( x[0] >= 2.20 ) && ( x[0] < 2.40 ) ) ans *= 0.980;
      // if ( ( x[0] >= 2.40 ) && ( x[0] < 2.60 ) ) ans *= 0.965;
      // if ( ( x[0] >= 2.60 ) && ( x[0] < 3.00 ) ) ans *= 0.950;
   } else if (fTrackingState == MpdPidUtils::kCF) {
      if ((x[0] > 0.10) && (x[0] < 0.15)) ans *= 1.044;
      if ((x[0] > 0.15) && (x[0] < 0.20)) ans *= 1.065;
      if ((x[0] > 0.20) && (x[0] < 0.25)) ans *= 1.072;
      if ((x[0] > 0.25) && (x[0] < 0.30)) ans *= 1.029;
      if ((x[0] > 0.30) && (x[0] < 0.45)) ans *= 0.990;
      if ((x[0] > 0.95) && (x[0] < 1.05)) ans *= 0.988;
      if ((x[0] > 1.10) && (x[0] < 1.60)) ans *= 1.009;
      if ((x[0] > 1.60) && (x[0] < 1.90)) ans *= 1.013;
      if ((x[0] > 1.90) && (x[0] < 2.20)) ans *= 1.021;
      if ((x[0] > 2.20) && (x[0] < 2.40)) ans *= 1.028;
      if ((x[0] > 2.40) && (x[0] < 2.60)) ans *= 1.035;
      if ((x[0] > 2.60) && (x[0] < 3.00)) ans *= 1.041;
   } else {
      if (x[0] < 0.15) ans *= 0.980;
      if (x[0] < 0.20) ans *= 0.970;
      if (x[0] < 0.25) ans *= 1.025;
      if (x[0] < 0.30) ans *= 1.020;
      if (x[0] < 0.40) ans *= 0.990;
      if (x[0] < 0.45) ans *= 0.990;
      if (x[0] < 0.50) ans *= 1.010;
      if (x[0] > 1.05) ans *= 1.015;
      if (x[0] > 1.50) ans *= 0.985;
      if (x[0] > 1.80) ans *= 0.990;
      if (x[0] > 2.10) ans *= 0.990;
      if (x[0] > 2.20) ans *= 0.990;
      if (x[0] > 2.40) ans *= 0.990;
      if (x[0] > 2.70) ans *= 0.985;
   }

   return ans;
}

Double_t MpdPid::parPrBB(Double_t *x, Double_t *par)
{
   Double_t x1, x2, x3, p[5], xint = fTrackingState != MpdPidUtils::kCFHM ? 0.384089 : 0.325918, ans;
   for (Int_t k = 0; k < 5; k++) p[k] = x[0] < xint ? par[k] : par[k + 5];
   x1  = p[0] / TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.88), p[3]);
   x2  = p[1] - TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.88), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / (x[0] / 0.9383), p[4]));
   ans = x1 * (x2 - x3);

   if (fTrackingState == MpdPidUtils::kCFHM) {
      if ((x[0] >= 0.05) && (x[0] < 0.10)) ans *= 0.86;
      if ((x[0] >= 0.10) && (x[0] < 0.15)) ans *= 0.975;
      if ((x[0] >= 0.15) && (x[0] < 0.20)) ans *= 1.04;
      if ((x[0] >= 0.20) && (x[0] < 0.25)) ans *= 1.015;
      if ((x[0] >= 0.30) && (x[0] < 0.40)) ans *= 0.98;
      // if ( ( x[0] >= 0.05 ) && ( x[0] < 0.10 ) ) ans *= 0.780;
      // if ( ( x[0] >= 0.10 ) && ( x[0] < 0.15 ) ) ans *= 0.925;
      // if ( ( x[0] >= 0.15 ) && ( x[0] < 0.20 ) ) ans *= 1.015;
      // if ( ( x[0] >= 0.20 ) && ( x[0] < 0.25 ) ) ans *= 0.985;
      // if ( ( x[0] >= 0.25 ) && ( x[0] < 0.30 ) ) ans *= 0.945;
      // if ( ( x[0] >= 0.30 ) && ( x[0] < 0.35 ) ) ans *= 0.970;
      // if ( ( x[0] >= 0.35 ) && ( x[0] < 0.40 ) ) ans *= 0.995;
      // if ( ( x[0] >= 2.20 ) && ( x[0] < 2.40 ) ) ans *= 0.995;
      // if ( ( x[0] >= 2.40 ) && ( x[0] < 2.60 ) ) ans *= 0.990;
      // if ( ( x[0] >= 2.60 ) && ( x[0] < 3.00 ) ) ans *= 0.975;
   } else if (fTrackingState == MpdPidUtils::kCF) {
      if ((x[0] > 0.05) && (x[0] < 0.10)) ans *= 0.891;
      if ((x[0] > 0.10) && (x[0] < 0.15)) ans *= 0.946;
      if ((x[0] > 0.15) && (x[0] < 0.20)) ans *= 1.050;
      if ((x[0] > 0.20) && (x[0] < 0.25)) ans *= 0.922;
      if ((x[0] > 0.25) && (x[0] < 0.30)) ans *= 0.848;
      if ((x[0] > 0.30) && (x[0] < 0.40)) ans *= 0.870;
      if ((x[0] > 0.45) && (x[0] < 0.55)) ans *= 0.988;
      if ((x[0] > 1.40) && (x[0] < 2.00)) ans *= 0.990;
      if ((x[0] > 2.10) && (x[0] < 2.20)) ans *= 1.010;
      if ((x[0] > 2.20) && (x[0] < 2.60)) ans *= 1.016;
      if ((x[0] > 2.60) && (x[0] < 3.00)) ans *= 1.024;
   } else {
      if (x[0] < 0.10) ans *= 0.690;
      if (x[0] < 0.15) ans *= 0.940;
      if (x[0] < 0.20) ans *= 0.955;
      if (x[0] < 0.25) ans *= 0.965;
      if (x[0] < 0.30) ans *= 0.975;
      if (x[0] < 0.45) ans *= 1.030;
      if (x[0] < 0.55) ans *= 0.970;
      if (x[0] < 0.60) ans *= 1.030;
      if (x[0] < 0.65) ans *= 0.975;
      if (x[0] < 0.70) ans *= 1.010;
      if (x[0] > 2.00) ans *= 1.010;
      if (x[0] > 2.40) ans *= 0.990;
   }
   return ans;
}

Double_t MpdPid::parDeBB(Double_t *x, Double_t *par)
{
   Double_t x1, x2, x3, p[5], xint = 0.40577186, ans;
   for (Int_t k = 0; k < 5; k++) p[k] = x[0] < xint ? par[k] : par[k + 5];
   x1  = p[0] / TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 3.52), p[3]);
   x2  = p[1] - TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 3.52), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / (x[0] / 1.876), p[4]));
   ans = x1 * (x2 - x3);

   if (fTrackingState == MpdPidUtils::kCFHM) {
      if ((x[0] >= 0.20) && (x[0] < 0.25)) ans *= 0.915;
      if ((x[0] >= 0.25) && (x[0] < 0.30)) ans *= 0.970;
      if ((x[0] >= 0.30) && (x[0] < 0.35)) ans *= 1.010;
      if ((x[0] >= 0.35) && (x[0] < 0.40)) ans *= 1.015;
      if ((x[0] >= 0.40) && (x[0] < 0.45)) ans *= 1.005;
      if ((x[0] >= 0.45) && (x[0] < 0.50)) ans *= 1.025;
      if ((x[0] >= 0.50) && (x[0] < 0.55)) ans *= 0.995;
      if ((x[0] >= 0.55) && (x[0] < 0.70)) ans *= 0.990;
   }
   return ans;
}

Double_t MpdPid::parTrBB(Double_t *x, Double_t *par)
{
   Double_t x1, x2, x3, p[5], xint = 0.53220409, ans;
   for (Int_t k = 0; k < 5; k++) p[k] = x[0] < xint ? par[k] : par[k + 5];
   x1  = p[0] / TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 7.89), p[3]);
   x2  = p[1] - TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 7.89), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / (x[0] / 2.81), p[4]));
   ans = x1 * (x2 - x3);

   if (fTrackingState == MpdPidUtils::kCFHM) {
      if ((x[0] >= 0.20) && (x[0] < 0.25)) ans *= 0.960;
      if ((x[0] >= 0.25) && (x[0] < 0.30)) ans *= 0.975;
      if ((x[0] >= 0.35) && (x[0] < 0.40)) ans *= 1.015;
      if ((x[0] >= 0.40) && (x[0] < 0.45)) ans *= 1.010;
      if ((x[0] >= 0.45) && (x[0] < 0.50)) ans *= 0.990;
      if ((x[0] >= 0.50) && (x[0] < 0.55)) ans *= 0.970;
      if ((x[0] >= 0.55) && (x[0] < 0.60)) ans *= 0.985;
      if ((x[0] >= 0.60) && (x[0] < 0.65)) ans *= 1.010;
      if ((x[0] >= 0.65) && (x[0] < 0.75)) ans *= 1.015;
      if ((x[0] >= 0.85) && (x[0] < 1.00)) ans *= 0.990;
   }
   return ans;
}

Double_t MpdPid::parHe3BB(Double_t *x, Double_t *par)
{
   Double_t x1, x2, x3, p[4], xint = 0.71388289, ans;
   for (Int_t k = 0; k < 4; k++) p[k] = x[0] < xint ? par[k] : par[k + 4];
   x1  = 1.0 + TMath::Power(x[0] / 1.4047, 2.0);
   x2  = TMath::Power(x[0] / 1.4047, p[3]);
   x3  = p[1] + p[2] * TMath::Log(1.0 + TMath::Power(x[0] / 1.4047, 2.0));
   ans = p[0] * (x1 / x2 * x3 - 1.0);

   if (fTrackingState == MpdPidUtils::kCFHM) {
      if ((x[0] >= 0.10) && (x[0] < 0.15)) ans *= 0.865;
      if ((x[0] >= 0.15) && (x[0] < 0.20)) ans *= 0.965;
      if ((x[0] >= 0.20) && (x[0] < 0.25)) ans *= 0.980;
      if ((x[0] >= 0.25) && (x[0] < 0.30)) ans *= 1.010;
      if ((x[0] >= 0.30) && (x[0] < 0.35)) ans *= 0.990;
      if ((x[0] >= 0.45) && (x[0] < 0.50)) ans *= 1.010;
      if ((x[0] >= 0.70) && (x[0] < 0.75)) ans *= 1.005;
      if ((x[0] >= 0.75) && (x[0] < 0.95)) ans *= 1.015;
      if ((x[0] >= 0.95) && (x[0] < 1.05)) ans *= 1.010;
      if ((x[0] >= 1.10) && (x[0] < 1.25)) ans *= 0.995;
      if ((x[0] >= 1.40) && (x[0] < 1.45)) ans *= 1.010;
      if ((x[0] >= 2.40) && (x[0] < 2.60)) ans *= 1.015;
   }
   return ans;
}

Double_t MpdPid::parHe4BB(Double_t *x, Double_t *par)
{
   Double_t x1, x2, x3, p[4], xint = 0.72405503, ans;
   for (Int_t k = 0; k < 4; k++) p[k] = x[0] < xint ? par[k] : par[k + 4];
   x1  = 1.0 + TMath::Power(x[0] / 1.863, 2.0);
   x2  = TMath::Power(x[0] / 1.863, p[3]);
   x3  = p[1] + p[2] * TMath::Log(1.0 + TMath::Power(x[0] / 1.863, 2.0));
   ans = p[0] * (x1 / x2 * x3 - 1.0);

   if (fTrackingState == MpdPidUtils::kCFHM) {
      if ((x[0] >= 0.10) && (x[0] < 0.15)) ans *= 0.755;
      if ((x[0] >= 0.15) && (x[0] < 0.20)) ans *= 0.850;
      if ((x[0] >= 0.20) && (x[0] < 0.25)) ans *= 0.950;
      if ((x[0] >= 0.30) && (x[0] < 0.35)) ans *= 1.015;
      if ((x[0] >= 0.65) && (x[0] < 0.70)) ans *= 1.010;
      if ((x[0] >= 0.70) && (x[0] < 0.75)) ans *= 1.015;
      if ((x[0] >= 0.75) && (x[0] < 0.80)) ans *= 1.040;
      if ((x[0] >= 0.80) && (x[0] < 0.85)) ans *= 1.065;
      if ((x[0] >= 0.85) && (x[0] < 1.05)) ans *= 1.080;
      if ((x[0] >= 1.05) && (x[0] < 1.10)) ans *= 1.065;
      if ((x[0] >= 1.10) && (x[0] < 1.20)) ans *= 1.050;
      if ((x[0] >= 1.20) && (x[0] < 1.25)) ans *= 1.040;
      if ((x[0] >= 1.25) && (x[0] < 1.30)) ans *= 1.030;
      if ((x[0] >= 1.30) && (x[0] < 1.40)) ans *= 1.015;
      if ((x[0] >= 1.60) && (x[0] < 2.20)) ans *= 0.995;
      if ((x[0] >= 2.40) && (x[0] < 2.60)) ans *= 1.030;
   }
   return ans;
}

Double_t MpdPid::AsymGaus(Double_t *x, Double_t *par)
{
   Double_t peak = par[1];
   if (x[0] < peak)
      return par[0] * TMath::Exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2));
   else
      return par[0] * TMath::Exp(-0.5 * TMath::Power(((x[0] - par[1]) / ((1. + par[3]) * par[2])), 2));
}

Double_t MpdPid::AsymGaus2(Double_t *x, Double_t *par)
{
   Double_t peak = par[1];
   if (x[0] < peak)
      return par[0] * TMath::Exp(-0.5 * TMath::Power(((x[0] - par[1]) / par[2]), 2)) *
             TMath::Exp(-0.5 * TMath::Power(((x[1] - par[4]) / par[5]), 2));
   else
      return par[0] * TMath::Exp(-0.5 * TMath::Power(((x[0] - par[1]) / ((1. + par[3]) * par[2])), 2)) *
             TMath::Exp(-0.5 * TMath::Power(((x[1] - par[4]) / par[5]), 2));
}

/// Sum of two exponents in total momentum (pions and electrons)
Double_t MpdPid::MomPi(Double_t *x, Double_t *par)
{
   Double_t p = x[0], xx, x1, x2;
   xx         = TMath::Sqrt(p * p + par[4] * par[4]) - par[4];
   x1         = (1.0 + par[1]) / (par[2] * (par[4] + par[2])) / TMath::Exp(xx / par[2]);
   x2         = par[1] / (par[3] * (par[4] + par[3])) / TMath::Exp(xx / par[3]);
   return (par[0] * p * (x1 + x2));
}

/// Difference of 2 exponents in total momentum (All species, except pi's and e+/-)
Double_t MpdPid::MomPr(Double_t *x, Double_t *par)
{
   Double_t p = x[0], xx, x1, x2;
   xx         = TMath::Sqrt(p * p + par[4] * par[4]) - par[4];
   x1         = (1.0 + par[1]) / (par[2] * (par[4] + par[2])) / TMath::Exp(xx / par[2]);
   x2         = par[1] / (par[3] * (par[4] + par[3])) / TMath::Exp(xx / par[3]);
   return (par[0] * p * (x1 - x2));
}

Double_t MpdPid::GetDedxWidthValue(Double_t p, MpdPidUtils::ePartType iType)
{
   Double_t WidthValue = 0.0;

   map<MpdPidUtils::ePartType, vecTF1ptrs>::iterator ret = fdEdxSigmaMap.find(iType);
   if (ret != fdEdxSigmaMap.end()) {
      for (TF1 *widthF : ret->second) {
         if ((p < widthF->GetXmax()) && (p >= widthF->GetXmin())) {
            WidthValue = widthF->Eval(p);
            break;
         }
      }
   }

   return WidthValue;
}

Double_t MpdPid::GetTailValue(Double_t p, MpdPidUtils::ePartType iType)
{
   Double_t TailValue = 0.0;

   map<MpdPidUtils::ePartType, vecTF1ptrs>::iterator ret = fdEdxDeltaMap.find(iType);
   if (ret != fdEdxDeltaMap.end()) {
      for (TF1 *asymF : ret->second) {
         if ((p < asymF->GetXmax()) && (p >= asymF->GetXmin())) {
            TailValue = asymF->Eval(p);
            break;
         }
      }
   }

   return TailValue;
}

Double_t MpdPid::ComputeEnLossSigma(Double_t p, Double_t dedx, Double_t emean, Double_t sige,
                                    MpdPidUtils::ePartType iType)
{
   Double_t sig             = TMath::Abs((dedx / emean - 1.) / sige);
   Double_t delta           = GetTailValue(p, iType);
   Double_t SideOfAsymmetry = delta * (dedx - emean);
   if (SideOfAsymmetry >= 0.0) sig /= 1. + TMath::Abs(delta);
   fEnLossSigmasArray[iType] = sig;
   return sig;
}

Double_t MpdPid::ComputeMSquaredSigma(Double_t m2, Double_t mmean, Double_t sigm, MpdPidUtils::ePartType iType)
{
   Double_t sig                = (m2 - mmean) / sigm;
   fMSquaredSigmasArray[iType] = TMath::Abs(sig);
   return sig;
}

Double_t MpdPid::ComputeDedxProb_asym(Double_t cut, Double_t p, Double_t dedx, Double_t n, Double_t emean,
                                      Double_t sige, MpdPidUtils::ePartType iType)
{
   Double_t Prob  = 0.0, xx;
   Double_t delta = GetTailValue(p, iType);
   xx             = ComputeEnLossSigma(p, dedx, emean, sige, iType);
   if (xx < cut) {
      fAsymGaus->SetParameters(n, 1., sige, delta);
      Prob = fAsymGaus->Eval(dedx / emean);
   }
   fMSquaredSigmasArray[iType] = -1.0;
   return Prob;
}

Double_t MpdPid::ComputeCombProb_asym(Double_t cut_dedx, Double_t cut_m2, Double_t p, Double_t dedx, Double_t m2,
                                      Double_t n, Double_t emean, Double_t mmean, Double_t sige, Double_t sigm,
                                      MpdPidUtils::ePartType iType)
{
   Double_t Prob  = 0.0, xx, yy;
   Double_t delta = GetTailValue(p, iType);
   xx             = ComputeEnLossSigma(p, dedx, emean, sige, iType);
   yy             = ComputeMSquaredSigma(m2, mmean, sigm, iType);
   if (iType == MpdPidUtils::kElectron) mmean = 0.002;
   if ((TMath::Abs(xx) < cut_dedx) && (TMath::Abs(yy) < cut_m2)) {
      fAsymGaus2->SetParameters(n, 1., sige, delta, mmean, sigm);
      Prob = fAsymGaus2->Eval(dedx / emean, m2);
   }
   return Prob;
}

void MpdPid::Init(TString Generator, TString Tracking, TString NSigPart, Double_t fCoef)
{
   cout << "MpdPid::Init().." << endl;

   Double_t PMIN = 0.0, PMAX = 5.0, xint;
   Double_t dedxParam;
   TString  sFunc;

   if (Tracking == "HP")
      fTrackingState = MpdPidUtils::kHP;
   else if (Tracking == "CF")
      fTrackingState = MpdPidUtils::kCF;
   else if (Tracking == "CFHM")
      fTrackingState = MpdPidUtils::kCFHM;
   else {
      cout << "ERROR! Unknown tracking method! Switch to default (\"CFHM\")." << endl;
      fTrackingState = MpdPidUtils::kCFHM;
   }

   /// Setting default ratio ('rat' is pos./neg.)
   fPrRatio = fEnergy < 7.0 ? 1000.0 : 100.0;

   vecTF1ptrs fBB[MpdPidUtils::kNSpecies], fm2[MpdPidUtils::kNSpecies], fSigDedx[MpdPidUtils::kNSpecies],
      fAsymDedx[MpdPidUtils::kNSpecies], fPartYield[MpdPidUtils::kNSpecies];

   ///
   /// Bethe-Bloch versus p
   ///

   TF1 *fParElBB = new TF1("fParElBB", this, &MpdPid::parElBB, PMIN, PMAX, 2, "MpdPid", "parElBB");
   fBB[MpdPidUtils::kElectron].push_back(fParElBB);
   fdEdxBBMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kElectron, fBB[MpdPidUtils::kElectron]));

   TF1 *fParMuBB = new TF1("fParMuBB", this, &MpdPid::parMuBB, PMIN, PMAX, 5, "MpdPid", "parMuBB");
   fBB[MpdPidUtils::kMuon].push_back(fParMuBB);
   fdEdxBBMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kMuon, fBB[MpdPidUtils::kMuon]));

   TF1 *fParPiBB = new TF1("fParPiBB", this, &MpdPid::parPiBB, PMIN, PMAX, 5, "MpdPid", "parPiBB");
   fBB[MpdPidUtils::kPion].push_back(fParPiBB);
   fdEdxBBMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fBB[MpdPidUtils::kPion]));

   TF1 *fParKaBB = new TF1("fParKaBB", this, &MpdPid::parKaBB, PMIN, PMAX, 10, "MpdPid", "parKaBB");
   fBB[MpdPidUtils::kKaon].push_back(fParKaBB);
   fdEdxBBMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fBB[MpdPidUtils::kKaon]));

   // xint = fTrackingState == MpdPidUtils::kCFHM ? 5.0 : 5.0;
   TF1 *fParPrBB = new TF1("fParPrBB", this, &MpdPid::parPrBB, PMIN, PMAX, 10, "MpdPid", "parPrBB");
   fBB[MpdPidUtils::kProton].push_back(fParPrBB);
   // if ( fTrackingState == MpdPidUtils::kCFHM ) {
   //	TF1* fParPrP1 = new TF1("fParPrP1","[0]",xint,5.0);
   //	fBB[MpdPidUtils::kProton].push_back(fParPrP1);
   // }
   fdEdxBBMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fBB[MpdPidUtils::kProton]));

   if (fTrackingState == MpdPidUtils::kCF) {
      xint          = 0.375;
      TF1 *fParDeP1 = new TF1("fParDeP1", "pol1(0)", PMIN, 0.225);
      fBB[MpdPidUtils::kDeuteron].push_back(fParDeP1);
      TF1 *fParDeP2 = new TF1("fParDeP2", "pol1(0)", 0.225, xint);
      fBB[MpdPidUtils::kDeuteron].push_back(fParDeP2);
   } else
      xint = PMIN;
   TF1 *fParDeBB = new TF1("fParDeBB", this, &MpdPid::parDeBB, xint, PMAX, 10, "MpdPid", "parDeBB");
   fBB[MpdPidUtils::kDeuteron].push_back(fParDeBB);
   fdEdxBBMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kDeuteron, fBB[MpdPidUtils::kDeuteron]));

   if (fTrackingState == MpdPidUtils::kCF) {
      xint          = 0.525;
      TF1 *fParTrP1 = new TF1("fParTrP1", "pol1(0)", PMIN, 0.3);
      fBB[MpdPidUtils::kTriton].push_back(fParTrP1);
      TF1 *fParTrP2 = new TF1("fParTrP2", "pol1(0)", 0.3, xint);
      fBB[MpdPidUtils::kTriton].push_back(fParTrP2);
   } else
      xint = PMIN;
   TF1 *fParTrBB = new TF1("fParTrBB", this, &MpdPid::parTrBB, xint, PMAX, 10, "MpdPid", "parTrBB");
   fBB[MpdPidUtils::kTriton].push_back(fParTrBB);
   fdEdxBBMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kTriton, fBB[MpdPidUtils::kTriton]));

   /// Double charged have p_reco= p_MC/2
   if (fTrackingState == MpdPidUtils::kCF) {
      xint           = 0.55;
      TF1 *fParHe3P1 = new TF1("fParHe3P1", "pol1(0)", PMIN, 0.24);
      fBB[MpdPidUtils::kHe3].push_back(fParHe3P1);
      TF1 *fParHe3P2 = new TF1("fParHe3P2", "pol1(0)", 0.24, 0.325);
      fBB[MpdPidUtils::kHe3].push_back(fParHe3P2);
      TF1 *fParHe3P3 = new TF1("fParHe3P3", "pol1(0)", 0.325, xint);
      fBB[MpdPidUtils::kHe3].push_back(fParHe3P3);
   } else
      xint = PMIN;
   PMAX           = fTrackingState == MpdPidUtils::kCFHM ? 2.6 : 5.0;
   TF1 *fParHe3BB = new TF1("fParHe3BB", this, &MpdPid::parHe3BB, xint, PMAX, 8, "MpdPid", "parHe3BB");
   fBB[MpdPidUtils::kHe3].push_back(fParHe3BB);
   PMAX = 5.0;
   if (fTrackingState == MpdPidUtils::kCFHM) {
      TF1 *fParHe3P4 = new TF1("fParHe3P4", "[0]", 2.6, PMAX);
      fBB[MpdPidUtils::kHe3].push_back(fParHe3P4);
   }
   fdEdxBBMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe3, fBB[MpdPidUtils::kHe3]));

   if (fTrackingState == MpdPidUtils::kCF) {
      xint           = 0.5;
      TF1 *fParHe4P1 = new TF1("fParHe4P1", "pol1(0)", PMIN, 0.32);
      fBB[MpdPidUtils::kHe4].push_back(fParHe4P1);
      TF1 *fParHe4P2 = new TF1("fParHe4P2", "pol1(0)", 0.32, xint);
      fBB[MpdPidUtils::kHe4].push_back(fParHe4P2);
   } else
      xint = PMIN;
   PMAX           = fTrackingState == MpdPidUtils::kCFHM ? 2.6 : 5.0;
   TF1 *fParHe4BB = new TF1("fParHe4BB", this, &MpdPid::parHe4BB, xint, PMAX, 8, "MpdPid", "parHe4BB");
   fBB[MpdPidUtils::kHe4].push_back(fParHe4BB);
   PMAX = 5.0;
   if (fTrackingState == MpdPidUtils::kCFHM) {
      TF1 *fParHe4P3 = new TF1("fParHe4P3", "[0]", 2.6, PMAX);
      fBB[MpdPidUtils::kHe4].push_back(fParHe4P3);
   }
   fdEdxBBMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe4, fBB[MpdPidUtils::kHe4]));

   ///
   /// Width of mass-squared versus total p
   ///

   TF1 *fParElM2 = new TF1("fParElM2", "pol2(0)", PMIN, PMAX);
   fm2[MpdPidUtils::kElectron].push_back(fParElM2);
   fParM2Map.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kElectron, fm2[MpdPidUtils::kElectron]));

   TF1 *fParMuM2 = new TF1("fParMuM2", "pol2(0)", PMIN, PMAX);
   fm2[MpdPidUtils::kMuon].push_back(fParMuM2);
   fParM2Map.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kMuon, fm2[MpdPidUtils::kMuon]));

   if (fTrackingState == MpdPidUtils::kCFHM) {
      TF1 *fParPiM2P1 = new TF1("fParPiM2P1", "[0]", PMIN, 0.5);
      fm2[MpdPidUtils::kPion].push_back(fParPiM2P1);
      PMIN = 0.5;
   }
   xint            = fTrackingState != MpdPidUtils::kCFHM ? 1.4 : 1.59086;
   TF1 *fParPiM2P2 = new TF1("fParPiM2P2", "pol2(0)", PMIN, xint);
   fm2[MpdPidUtils::kPion].push_back(fParPiM2P2);
   TF1 *fParPiM2P3 = new TF1("fParPiM2P3", "pol2(0)", xint, PMAX);
   fm2[MpdPidUtils::kPion].push_back(fParPiM2P3);
   fParM2Map.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fm2[MpdPidUtils::kPion]));
   PMIN = 0.0;

   if (fTrackingState != MpdPidUtils::kCFHM) {
      TF1 *fParKaM2P1 = new TF1("fParKaM2P1", "pol2(0)", PMIN, PMAX);
      fm2[MpdPidUtils::kKaon].push_back(fParKaM2P1);
   } else {
      TF1 *fParKaM2P1 = new TF1("fParKaM2P1", "pol3(0)", PMIN, PMAX);
      fm2[MpdPidUtils::kKaon].push_back(fParKaM2P1);
      // TF1* fParKaM2P1 = new TF1("fParKaM2P1","[0]",PMIN,0.2);
      // fm2[MpdPidUtils::kKaon].push_back(fParKaM2P1);
      // TF1* fParKaM2P2 = new TF1("fParKaM2P2","pol2(0)",0.2,0.705661);
      // fm2[MpdPidUtils::kKaon].push_back(fParKaM2P2);
      // TF1* fParKaM2P3 = new TF1("fParKaM2P3","pol2(0)",0.705661,PMAX);
      // fm2[MpdPidUtils::kKaon].push_back(fParKaM2P3);
   }
   fParM2Map.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fm2[MpdPidUtils::kKaon]));

   if (fTrackingState == MpdPidUtils::kCFHM) {
      TF1 *fParPrM2P1 = new TF1("fParPrM2P1", "[0]", PMIN, 0.3);
      fm2[MpdPidUtils::kProton].push_back(fParPrM2P1);
      PMIN = 0.3;
   }
   xint            = fTrackingState != MpdPidUtils::kCFHM ? 1.4 : 0.930546;
   sFunc           = fTrackingState != MpdPidUtils::kCFHM ? "pol3(0)" : "pol2(0)";
   TF1 *fParPrM2P2 = new TF1("fParPrM2P2", sFunc, PMIN, xint);
   fm2[MpdPidUtils::kProton].push_back(fParPrM2P2);
   TF1 *fParPrM2P3 = new TF1("fParPrM2P3", "pol2(0)", xint, PMAX);
   fm2[MpdPidUtils::kProton].push_back(fParPrM2P3);
   fParM2Map.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fm2[MpdPidUtils::kProton]));
   PMIN = 0.0;

   if (fTrackingState == MpdPidUtils::kCFHM) {
      TF1 *fParDeM2P1 = new TF1("fParDeM2P1", "[0]", PMIN, 0.5);
      fm2[MpdPidUtils::kDeuteron].push_back(fParDeM2P1);
      PMIN = 0.5;
   }
   sFunc         = fTrackingState != MpdPidUtils::kCFHM ? "pol3(0)" : "pol2(0)";
   TF1 *fParDeM2 = new TF1("fParDeM2", sFunc, PMIN, PMAX);
   fm2[MpdPidUtils::kDeuteron].push_back(fParDeM2);
   fParM2Map.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kDeuteron, fm2[MpdPidUtils::kDeuteron]));
   PMIN = 0.0;

   if (fTrackingState == MpdPidUtils::kCFHM) {
      TF1 *fParTrM2P1 = new TF1("fParTrM2P1", "[0]", PMIN, 0.7);
      fm2[MpdPidUtils::kTriton].push_back(fParTrM2P1);
      PMIN = 0.7;
   }
   sFunc         = fTrackingState == MpdPidUtils::kHP ? "pol3(0)" : "pol2(0)";
   TF1 *fParTrM2 = new TF1("fParTrM2", sFunc, PMIN, PMAX);
   fm2[MpdPidUtils::kTriton].push_back(fParTrM2);
   fParM2Map.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kTriton, fm2[MpdPidUtils::kTriton]));
   PMIN = 0.0;

   if (fTrackingState == MpdPidUtils::kCFHM) {
      TF1 *fParHe3M2P1 = new TF1("fParHe3M2P1", "[0]", PMIN, 0.55);
      fm2[MpdPidUtils::kHe3].push_back(fParHe3M2P1);
      PMIN = 0.55;
   }
   sFunc          = fTrackingState != MpdPidUtils::kCFHM ? "[0]" : "pol2(0)";
   TF1 *fParHe3M2 = new TF1("fParHe3M2", sFunc, PMIN, PMAX);
   fm2[MpdPidUtils::kHe3].push_back(fParHe3M2);
   fParM2Map.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe3, fm2[MpdPidUtils::kHe3]));
   PMIN = 0.0;

   if (fTrackingState == MpdPidUtils::kCFHM) {
      TF1 *fParHe4M2P1 = new TF1("fParHe4M2P1", "[0]", PMIN, 0.65);
      fm2[MpdPidUtils::kHe4].push_back(fParHe4M2P1);
      PMIN = 0.65;
   }
   sFunc          = fTrackingState != MpdPidUtils::kCFHM ? "[0]" : "pol2(0)";
   TF1 *fParHe4M2 = new TF1("fParHe4M2", sFunc, PMIN, PMAX);
   fm2[MpdPidUtils::kHe4].push_back(fParHe4M2);
   fParM2Map.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe4, fm2[MpdPidUtils::kHe4]));
   PMIN = 0.0;

   map<MpdPidUtils::ePartType, vecTF1ptrs>::iterator it;
   if (fTrackingState == MpdPidUtils::kHP) {
      /// dE/dx Bethe-Bloch
      it = fdEdxBBMap.find(MpdPidUtils::kElectron);
      it->second[0]->SetParameters(fCoef * (-8.27289e-09), fCoef * 2.10438e-06);
      it = fdEdxBBMap.find(MpdPidUtils::kMuon);
      it->second[0]->SetParameters(fCoef * 2.56183e-06, 2.7034, 2.52734, 0.854563, 0.330576);
      it = fdEdxBBMap.find(MpdPidUtils::kPion);
      it->second[0]->SetParameters(fCoef * (-1.19342e-07), -7.83114, 8.17749, 1.85775, -1.80695);
      it = fdEdxBBMap.find(MpdPidUtils::kKaon);
      it->second[0]->SetParameters(fCoef * 6.50167e-07, 1.01718, -0.795357, 1.80916, 0.0707667, fCoef * 6.50167e-07,
                                   1.01718, -0.795357, 1.80916, 0.0707667);
      it = fdEdxBBMap.find(MpdPidUtils::kProton);
      it->second[0]->SetParameters(fCoef * 4.40008e-07, 2.97563, -0.192657, 2.16118, 0.61451, fCoef * 4.40008e-07,
                                   2.97563, -0.192657, 2.16118, 0.61451);
      it = fdEdxBBMap.find(MpdPidUtils::kDeuteron);
      it->second[0]->SetParameters(fCoef * 3.27e-07, 3.74, -0.23, 2.32, 0.987, fCoef * 3.27e-07, 3.74, -0.23, 2.32,
                                   0.987);
      it = fdEdxBBMap.find(MpdPidUtils::kTriton);
      it->second[0]->SetParameters(fCoef * 2.59e-07, 5.06, 0.0001, 2.2, 1.056, fCoef * 2.59e-07, 5.06, 0.0001, 2.2,
                                   1.056);
      it = fdEdxBBMap.find(MpdPidUtils::kHe3);
      it->second[0]->SetParameters(fCoef * 2.86201e-06, 2.10168, 2.74807e-01, 1.86774, fCoef * 2.86201e-06, 2.10168,
                                   2.74807e-01, 1.86774);
      it = fdEdxBBMap.find(MpdPidUtils::kHe4);
      it->second[0]->SetParameters(fCoef * 2.96e-06, 2.085, 0.256, 1.85, fCoef * 2.96e-06, 2.085, 0.256, 1.85);

      /// Width of m^2
      it = fParM2Map.find(MpdPidUtils::kElectron);
      it->second[0]->SetParameters(0.00102552, -0.000243946, 0.0307395);
      it = fParM2Map.find(MpdPidUtils::kMuon);
      it->second[0]->SetParameters(0.00155957, -0.000984273, 0.0306857);
      it = fParM2Map.find(MpdPidUtils::kPion);
      it->second[0]->SetParameters(0.00259115, -0.00251021, 0.0287287);
      it->second[1]->SetParameters(-0.0393955, 0.0586315, 0.00516675);
      it = fParM2Map.find(MpdPidUtils::kKaon);
      it->second[0]->SetParameters(0.00144014, 0.0183536, 0.0161613);
      it = fParM2Map.find(MpdPidUtils::kProton);
      it->second[0]->SetParameters(0.0777042, -0.123828, 0.139278, -0.0338542);
      it->second[1]->SetParameters(0.0244298, 0.018534, 0.0174998);
      it = fParM2Map.find(MpdPidUtils::kDeuteron);
      it->second[0]->SetParameters(0.535691, -0.529882, 0.293807, -0.0428139);
      it = fParM2Map.find(MpdPidUtils::kTriton);
      it->second[0]->SetParameters(0.422, 0.3, -0.202, 0.0524);
      it = fParM2Map.find(MpdPidUtils::kHe3);
      it->second[0]->SetParameter(0, 0.17);
      it = fParM2Map.find(MpdPidUtils::kHe4);
      it->second[0]->SetParameter(0, 0.3);

      const Double_t kSigmaDedx[3] = {0.05, 0.11, 0.08}; ///< constant values of dE/dx width (Hit Producer Tracking)
      /// Width of dE/dx
      TF1 *fEnLossSigmaEl = new TF1("fEnLossSigmaEl", "[0]", PMIN, PMAX);
      fEnLossSigmaEl->SetParameter(0, kSigmaDedx[2]);
      fSigDedx[MpdPidUtils::kElectron].push_back(fEnLossSigmaEl);
      fdEdxSigmaMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kElectron, fSigDedx[MpdPidUtils::kElectron]));

      TF1 *fEnLossSigmaMuP1 = new TF1("fEnLossSigmaMuP1", "[0]", PMIN, 0.15);
      fEnLossSigmaMuP1->SetParameter(0, kSigmaDedx[1]);
      TF1 *fEnLossSigmaMuP2 = new TF1("fEnLossSigmaMuP2", "[0]", 0.15, PMAX);
      fEnLossSigmaMuP2->SetParameter(0, kSigmaDedx[0]);
      fSigDedx[MpdPidUtils::kMuon].push_back(fEnLossSigmaMuP1);
      fSigDedx[MpdPidUtils::kMuon].push_back(fEnLossSigmaMuP2);
      fdEdxSigmaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kMuon, fSigDedx[MpdPidUtils::kMuon]));

      TF1 *fEnLossSigmaPiP1 = new TF1("fEnLossSigmaPiP1", "[0]", PMIN, 0.15);
      fEnLossSigmaPiP1->SetParameter(0, kSigmaDedx[1]);
      TF1 *fEnLossSigmaPiP2 = new TF1("fEnLossSigmaPiP2", "[0]", 0.15, PMAX);
      fEnLossSigmaPiP2->SetParameter(0, kSigmaDedx[0]);
      fSigDedx[MpdPidUtils::kPion].push_back(fEnLossSigmaPiP1);
      fSigDedx[MpdPidUtils::kPion].push_back(fEnLossSigmaPiP2);
      fdEdxSigmaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fSigDedx[MpdPidUtils::kPion]));

      TF1 *fEnLossSigmaKaP1 = new TF1("fEnLossSigmaKaP1", "[0]", PMIN, 0.15);
      fEnLossSigmaKaP1->SetParameter(0, kSigmaDedx[1]);
      TF1 *fEnLossSigmaKaP2 = new TF1("fEnLossSigmaKaP2", "[0]", 0.15, PMAX);
      fEnLossSigmaKaP2->SetParameter(0, kSigmaDedx[0]);
      fSigDedx[MpdPidUtils::kKaon].push_back(fEnLossSigmaKaP1);
      fSigDedx[MpdPidUtils::kKaon].push_back(fEnLossSigmaKaP2);
      fdEdxSigmaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fSigDedx[MpdPidUtils::kKaon]));

      TF1 *fEnLossSigmaPrP1 = new TF1("fEnLossSigmaPrP1", "[0]", PMIN, 0.15);
      fEnLossSigmaPrP1->SetParameter(0, kSigmaDedx[1]);
      TF1 *fEnLossSigmaPrP2 = new TF1("fEnLossSigmaPrP2", "[0]", 0.15, PMAX);
      fEnLossSigmaPrP2->SetParameter(0, kSigmaDedx[0]);
      fSigDedx[MpdPidUtils::kProton].push_back(fEnLossSigmaPrP1);
      fSigDedx[MpdPidUtils::kProton].push_back(fEnLossSigmaPrP2);
      fdEdxSigmaMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fSigDedx[MpdPidUtils::kProton]));

      TF1 *fEnLossSigmaDeP1 = new TF1("fEnLossSigmaDeP1", "[0]", PMIN, 0.45);
      fEnLossSigmaDeP1->SetParameter(0, kSigmaDedx[1]);
      TF1 *fEnLossSigmaDeP2 = new TF1("fEnLossSigmaDeP2", "[0]", 0.45, PMAX);
      fEnLossSigmaDeP2->SetParameter(0, kSigmaDedx[0]);
      fSigDedx[MpdPidUtils::kDeuteron].push_back(fEnLossSigmaDeP1);
      fSigDedx[MpdPidUtils::kDeuteron].push_back(fEnLossSigmaDeP2);
      fdEdxSigmaMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kDeuteron, fSigDedx[MpdPidUtils::kDeuteron]));

      TF1 *fEnLossSigmaTrP1 = new TF1("fEnLossSigmaTrP1", "[0]", PMIN, 0.45);
      fEnLossSigmaTrP1->SetParameter(0, kSigmaDedx[1]);
      TF1 *fEnLossSigmaTrP2 = new TF1("fEnLossSigmaTrP2", "[0]", 0.45, PMAX);
      fEnLossSigmaTrP2->SetParameter(0, kSigmaDedx[0]);
      fSigDedx[MpdPidUtils::kTriton].push_back(fEnLossSigmaTrP1);
      fSigDedx[MpdPidUtils::kTriton].push_back(fEnLossSigmaTrP2);
      fdEdxSigmaMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kTriton, fSigDedx[MpdPidUtils::kTriton]));

      TF1 *fEnLossSigmaHe3P1 = new TF1("fEnLossSigmaHe3P1", "[0]", PMIN, 0.45);
      fEnLossSigmaHe3P1->SetParameter(0, kSigmaDedx[1]);
      TF1 *fEnLossSigmaHe3P2 = new TF1("fEnLossSigmaHe3P2", "[0]", 0.45, PMAX);
      fEnLossSigmaHe3P2->SetParameter(0, kSigmaDedx[0]);
      fSigDedx[MpdPidUtils::kHe3].push_back(fEnLossSigmaHe3P1);
      fSigDedx[MpdPidUtils::kHe3].push_back(fEnLossSigmaHe3P2);
      fdEdxSigmaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe3, fSigDedx[MpdPidUtils::kHe3]));

      TF1 *fEnLossSigmaHe4P1 = new TF1("fEnLossSigmaHe4P1", "[0]", PMIN, PMAX);
      fEnLossSigmaHe4P1->SetParameter(0, kSigmaDedx[1]);
      TF1 *fEnLossSigmaHe4P2 = new TF1("fEnLossSigmaHe4P2", "[0]", PMIN, PMAX);
      fEnLossSigmaHe4P2->SetParameter(0, kSigmaDedx[0]);
      fSigDedx[MpdPidUtils::kHe4].push_back(fEnLossSigmaHe4P1);
      fSigDedx[MpdPidUtils::kHe4].push_back(fEnLossSigmaHe4P2);
      fdEdxSigmaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe4, fSigDedx[MpdPidUtils::kHe4]));

      /// Asymmetry parameter of dE/dx
      TF1 *fAsymElP1 = new TF1("fAsymElP1", "pol1(0)", PMIN, 0.5);
      fAsymElP1->SetParameters(0.456871, -0.0424708);
      TF1 *fAsymElP2 = new TF1("fAsymElP2", "pol1(0)", 0.5, 0.85);
      fAsymElP2->SetParameters(0.4418036, -0.0131662);
      TF1 *fAsymElP3 = new TF1("fAsymElP3", "pol1(0)", 0.85, 3.0);
      fAsymElP3->SetParameters(0.440012, -0.0113746);
      TF1 *fAsymElP4 = new TF1("fAsymElP4", "[0]", 3.0, PMAX);
      fAsymElP4->SetParameter(0, fAsymElP3->Eval(3.0));
      fAsymDedx[MpdPidUtils::kElectron].push_back(fAsymElP1);
      fAsymDedx[MpdPidUtils::kElectron].push_back(fAsymElP2);
      fAsymDedx[MpdPidUtils::kElectron].push_back(fAsymElP3);
      fAsymDedx[MpdPidUtils::kElectron].push_back(fAsymElP4);
      fdEdxDeltaMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kElectron, fAsymDedx[MpdPidUtils::kElectron]));

      TF1 *fAsymMuP1 = new TF1("fAsymMuP1", "pol1(0)", PMIN, 0.175);
      fAsymMuP1->SetParameters(12.138, -66.6292);
      TF1 *fAsymMuP2 = new TF1("fAsymMuP2", "pol1(0)", 0.175, 0.65);
      fAsymMuP2->SetParameters(0.31412, 0.158834);
      TF1 *fAsymMuP3 = new TF1("fAsymMuP3", "pol1(0)", 0.65, 3.0);
      fAsymMuP3->SetParameters(0.409629, -0.0234975);
      TF1 *fAsymMuP4 = new TF1("fAsymMuP4", "[0]", 3.0, PMAX);
      fAsymMuP4->SetParameter(0, fAsymMuP3->Eval(3.0));
      fAsymDedx[MpdPidUtils::kMuon].push_back(fAsymMuP1);
      fAsymDedx[MpdPidUtils::kMuon].push_back(fAsymMuP2);
      fAsymDedx[MpdPidUtils::kMuon].push_back(fAsymMuP3);
      fAsymDedx[MpdPidUtils::kMuon].push_back(fAsymMuP4);
      fdEdxDeltaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kMuon, fAsymDedx[MpdPidUtils::kMuon]));

      TF1 *fAsymPiP1 = new TF1("fAsymPiP1", "pol1(0)", PMIN, 0.18);
      fAsymPiP1->SetParameters(11.3953, -58.1321);
      TF1 *fAsymPiP2 = new TF1("fAsymPiP2", "pol1(0)", 0.18, 0.9);
      fAsymPiP2->SetParameters(0.412467, -0.0289088);
      TF1 *fAsymPiP3 = new TF1("fAsymPiP3", "pol1(0)", 0.9, 3.0);
      fAsymPiP3->SetParameters(0.408242, -0.0520603);
      TF1 *fAsymPiP4 = new TF1("fAsymPiP4", "[0]", 3.0, PMAX);
      fAsymPiP4->SetParameter(0, fAsymPiP3->Eval(3.0));
      fAsymDedx[MpdPidUtils::kPion].push_back(fAsymPiP1);
      fAsymDedx[MpdPidUtils::kPion].push_back(fAsymPiP2);
      fAsymDedx[MpdPidUtils::kPion].push_back(fAsymPiP3);
      fAsymDedx[MpdPidUtils::kPion].push_back(fAsymPiP4);
      fdEdxDeltaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fAsymDedx[MpdPidUtils::kPion]));

      TF1 *fAsymKaP1 = new TF1("fAsymKaP1", "pol1(0)", PMIN, 0.45);
      fAsymKaP1->SetParameters(3.21698, -6.37965);
      TF1 *fAsymKaP2 = new TF1("fAsymKaP2", "pol1(0)", 0.45, 1.05);
      fAsymKaP2->SetParameters(0.237053, 0.143231);
      TF1 *fAsymKaP3 = new TF1("fAsymKaP3", "pol1(0)", 1.05, 3.0);
      fAsymKaP3->SetParameters(0.252239, 0.0558685);
      TF1 *fAsymKaP4 = new TF1("fAsymKaP4", "[0]", 3.0, PMAX);
      fAsymKaP4->SetParameter(0, fAsymKaP3->Eval(3.0));
      fAsymDedx[MpdPidUtils::kKaon].push_back(fAsymKaP1);
      fAsymDedx[MpdPidUtils::kKaon].push_back(fAsymKaP2);
      fAsymDedx[MpdPidUtils::kKaon].push_back(fAsymKaP3);
      fAsymDedx[MpdPidUtils::kKaon].push_back(fAsymKaP4);
      fdEdxDeltaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fAsymDedx[MpdPidUtils::kKaon]));

      TF1 *fAsymPrP1 = new TF1("fAsymPrP1", "pol1(0)", PMIN, 0.3);
      fAsymPrP1->SetParameters(6.40228, -19.5847);
      TF1 *fAsymPrP2 = new TF1("fAsymPrP2", "pol1(0)", 0.3, 1.0);
      fAsymPrP2->SetParameters(0.193311, 0.179621);
      TF1 *fAsymPrP3 = new TF1("fAsymPrP3", "pol1(0)", 1.0, 3.0);
      fAsymPrP3->SetParameters(0.21947, 0.0936557);
      TF1 *fAsymPrP4 = new TF1("fAsymPrP4", "[0]", 3.0, PMAX);
      fAsymPrP4->SetParameter(0, fAsymPrP3->Eval(3.0));
      fAsymDedx[MpdPidUtils::kProton].push_back(fAsymPrP1);
      fAsymDedx[MpdPidUtils::kProton].push_back(fAsymPrP2);
      fAsymDedx[MpdPidUtils::kProton].push_back(fAsymPrP3);
      fAsymDedx[MpdPidUtils::kProton].push_back(fAsymPrP4);
      fdEdxDeltaMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fAsymDedx[MpdPidUtils::kProton]));
   } else if (fTrackingState == MpdPidUtils::kCF) {
      /// dE/dx Bethe-Bloch
      it = fdEdxBBMap.find(MpdPidUtils::kElectron);
      it->second[0]->SetParameters(fCoef * (-65.7432), fCoef * 4007.68);
      it = fdEdxBBMap.find(MpdPidUtils::kMuon);
      it->second[0]->SetParameters(fCoef * (-54.6053), -41.38, 109.594, 1.30576, -4.66578);
      it = fdEdxBBMap.find(MpdPidUtils::kPion);
      it->second[0]->SetParameters(fCoef * (-345.254), -4.65898, 3.41161, 1.88069, -1.20667);
      it = fdEdxBBMap.find(MpdPidUtils::kKaon);
      it->second[0]->SetParameters(fCoef * (-665.977), -4.10879, 0.393422, 1.01384, 4.12062, fCoef * (-665.977),
                                   -4.10879, 0.393422, 1.01384, 4.12062);
      it = fdEdxBBMap.find(MpdPidUtils::kProton);
      it->second[0]->SetParameters(fCoef * (-6140.11), 0.8323, 2.03597, 1.34309, 0.746625, fCoef * (-16259.2), 3.08434,
                                   9.53496, -0.128409, 2.79948);
      it = fdEdxBBMap.find(MpdPidUtils::kDeuteron);
      it->second[0]->SetParameters(fCoef * (112.e+04), fCoef * (-4466.666e+03));
      it->second[1]->SetParameters(fCoef * (215507.), fCoef * (-474223.));
      it->second[2]->SetParameters(fCoef * (-3591.57), 0.656424, 1.61966, 0.585121, 2.70803, fCoef * (-3591.57),
                                   0.656424, 1.61966, 0.585121, 2.70803);
      it = fdEdxBBMap.find(MpdPidUtils::kTriton);
      it->second[0]->SetParameters(fCoef * (1.094542e+06), fCoef * (-1.81963e+06));
      it->second[1]->SetParameters(fCoef * (226703.), fCoef * (-362129.));
      it->second[2]->SetParameters(fCoef * (-6467.59), 1.83428, 3.94612, 0.0819389, 3.2898, fCoef * (-6467.59), 1.83428,
                                   3.94612, 0.0819389, 3.2898);
      it = fdEdxBBMap.find(MpdPidUtils::kHe3);
      it->second[0]->SetParameters(fCoef * (3.09e+06), fCoef * (-10.833333e+06));
      it->second[1]->SetParameters(fCoef * (1.403646e+06), fCoef * (-3.806859e+06));
      it->second[2]->SetParameters(fCoef * (333495.), fCoef * (-514085.));
      it->second[3]->SetParameters(fCoef * (-19120.8), -0.183431, 0.323092, 2.41984, fCoef * (-19120.8), -0.183431,
                                   0.323092, 2.41984);
      it = fdEdxBBMap.find(MpdPidUtils::kHe4);
      it->second[0]->SetParameters(fCoef * (3.243636e+06), fCoef * (-8.636364e+06));
      it->second[1]->SetParameters(fCoef * (1.100082e+06), fCoef * (-1.937756e+06));
      it->second[2]->SetParameters(fCoef * (-16150.2), -0.0645649, -0.0213786, 3.61265, fCoef * (-16150.2), -0.0645649,
                                   -0.0213786, 3.61265);

      /// m2
      it = fParM2Map.find(MpdPidUtils::kElectron);
      it->second[0]->SetParameters(0.001227, -0.000973509, 0.0314155);
      it = fParM2Map.find(MpdPidUtils::kMuon);
      it->second[0]->SetParameters(0.00166279, -0.00131341, 0.0311028);
      it = fParM2Map.find(MpdPidUtils::kPion);
      it->second[0]->SetParameters(0.00284751, -0.00311967, 0.0238384);
      it->second[1]->SetParameters(-0.0139431, 0.0254295, 0.0111476);
      it = fParM2Map.find(MpdPidUtils::kKaon);
      it->second[0]->SetParameters(0.00520284, 0.0092726, 0.0177252);
      it = fParM2Map.find(MpdPidUtils::kProton);
      it->second[0]->SetParameters(0.0873237, -0.127771, 0.130356, -0.0290419);
      it->second[1]->SetParameters(0.00606334, 0.040014, 0.0113814);
      it = fParM2Map.find(MpdPidUtils::kDeuteron);
      it->second[0]->SetParameters(0.417362, -0.399356, 0.24154, -0.0339673);
      it = fParM2Map.find(MpdPidUtils::kTriton);
      it->second[0]->SetParameters(1.03599, -0.663541, 0.199024);
      it = fParM2Map.find(MpdPidUtils::kHe3);
      it->second[0]->SetParameter(0, 0.13);
      it = fParM2Map.find(MpdPidUtils::kHe4);
      it->second[0]->SetParameter(0, 0.24);

      /// dE/dx width
      TF1 *fEnLossSigmaElP1 = new TF1("fEnLossSigmaElP1", "pol1(0)", PMIN, 0.55);
      fEnLossSigmaElP1->SetParameters(0.0747217, -0.0308101);
      TF1 *fEnLossSigmaElP2 = new TF1("fEnLossSigmaElP2", "pol1(0)", 0.55, 2.0);
      fEnLossSigmaElP2->SetParameters(0.0653074, -0.00546004);
      TF1 *fEnLossSigmaElP3 = new TF1("fEnLossSigmaElP3", "pol1(0)", 2.0, 3.0);
      fEnLossSigmaElP3->SetParameters(0.0572145, -0.00104922);
      TF1 *fEnLossSigmaElP4 = new TF1("fEnLossSigmaElP4", "[0]", 3.0, PMAX);
      fEnLossSigmaElP4->SetParameter(0, fEnLossSigmaElP3->Eval(3.0));
      fSigDedx[MpdPidUtils::kElectron].push_back(fEnLossSigmaElP1);
      fSigDedx[MpdPidUtils::kElectron].push_back(fEnLossSigmaElP2);
      fSigDedx[MpdPidUtils::kElectron].push_back(fEnLossSigmaElP3);
      fSigDedx[MpdPidUtils::kElectron].push_back(fEnLossSigmaElP4);
      fdEdxSigmaMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kElectron, fSigDedx[MpdPidUtils::kElectron]));

      TF1 *fEnLossSigmaMuP1 = new TF1("fEnLossSigmaMuP1", "pol1(0)", PMIN, 0.3);
      fEnLossSigmaMuP1->SetParameters(0.115611, -0.179933);
      TF1 *fEnLossSigmaMuP2 = new TF1("fEnLossSigmaMuP2", "pol1(0)", 0.3, 1.05);
      fEnLossSigmaMuP2->SetParameters(0.0583018, 0.00140717);
      TF1 *fEnLossSigmaMuP3 = new TF1("fEnLossSigmaMuP3", "pol1(0)", 1.05, 3.0);
      fEnLossSigmaMuP3->SetParameters(0.0570477, -8.02801e-05);
      TF1 *fEnLossSigmaMuP4 = new TF1("fEnLossSigmaMuP4", "[0]", 3.0, PMAX);
      fEnLossSigmaMuP4->SetParameter(0, fEnLossSigmaMuP3->Eval(3.0));
      fSigDedx[MpdPidUtils::kMuon].push_back(fEnLossSigmaMuP1);
      fSigDedx[MpdPidUtils::kMuon].push_back(fEnLossSigmaMuP2);
      fSigDedx[MpdPidUtils::kMuon].push_back(fEnLossSigmaMuP3);
      fSigDedx[MpdPidUtils::kMuon].push_back(fEnLossSigmaMuP4);
      fdEdxSigmaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kMuon, fSigDedx[MpdPidUtils::kMuon]));

      TF1 *fEnLossSigmaPiP1 = new TF1("fEnLossSigmaPiP1", "pol1(0)", PMIN, 0.319058);
      fEnLossSigmaPiP1->SetParameters(0.103709, -0.142807);
      TF1 *fEnLossSigmaPiP2 = new TF1("fEnLossSigmaPiP2", "pol1(0)", 0.319058, 1.41097);
      fEnLossSigmaPiP2->SetParameters(0.0568649, 0.0040131);
      TF1 *fEnLossSigmaPiP3 = new TF1("fEnLossSigmaPiP3", "pol1(0)", 1.41097, 1.6955);
      fEnLossSigmaPiP3->SetParameters(0.0606278, 0.00134621);
      TF1 *fEnLossSigmaPiP4 = new TF1("fEnLossSigmaPiP4", "pol1(0)", 1.6855, 3.0);
      fEnLossSigmaPiP4->SetParameters(0.0527141, 0.00604137);
      TF1 *fEnLossSigmaPiP5 = new TF1("fEnLossSigmaPiP5", "[0]", 3.0, PMAX);
      fEnLossSigmaPiP5->SetParameter(0, fEnLossSigmaPiP4->Eval(3.0));
      fSigDedx[MpdPidUtils::kPion].push_back(fEnLossSigmaPiP1);
      fSigDedx[MpdPidUtils::kPion].push_back(fEnLossSigmaPiP2);
      fSigDedx[MpdPidUtils::kPion].push_back(fEnLossSigmaPiP3);
      fSigDedx[MpdPidUtils::kPion].push_back(fEnLossSigmaPiP4);
      fSigDedx[MpdPidUtils::kPion].push_back(fEnLossSigmaPiP5);
      fdEdxSigmaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fSigDedx[MpdPidUtils::kPion]));

      TF1 *fEnLossSigmaKaP1 = new TF1("fEnLossSigmaKaP1", "pol1(0)", PMIN, 0.202208);
      fEnLossSigmaKaP1->SetParameters(0.122114, -0.299878);
      TF1 *fEnLossSigmaKaP2 = new TF1("fEnLossSigmaKaP2", "pol1(0)", 0.202208, 3.0);
      fEnLossSigmaKaP2->SetParameters(0.0614433, 0.000163801);
      TF1 *fEnLossSigmaKaP3 = new TF1("fEnLossSigmaKaP3", "[0]", 3.0, PMAX);
      fEnLossSigmaKaP3->SetParameter(0, fEnLossSigmaKaP2->Eval(3.0));
      fSigDedx[MpdPidUtils::kKaon].push_back(fEnLossSigmaKaP1);
      fSigDedx[MpdPidUtils::kKaon].push_back(fEnLossSigmaKaP2);
      fSigDedx[MpdPidUtils::kKaon].push_back(fEnLossSigmaKaP3);
      fdEdxSigmaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fSigDedx[MpdPidUtils::kKaon]));

      TF1 *fEnLossSigmaPrP1 = new TF1("fEnLossSigmaPrP1", "pol1(0)", PMIN, 0.5);
      fEnLossSigmaPrP1->SetParameters(0.134444, -0.148889);
      TF1 *fEnLossSigmaPrP2 = new TF1("fEnLossSigmaPrP2", "pol1(0)", 0.5, 3.0);
      fEnLossSigmaPrP2->SetParameters(0.0623534, -0.00287991);
      TF1 *fEnLossSigmaPrP3 = new TF1("fEnLossSigmaPrP3", "[0]", 3.0, PMAX);
      fEnLossSigmaPrP3->SetParameter(0, fEnLossSigmaPrP2->Eval(3.0));
      fSigDedx[MpdPidUtils::kProton].push_back(fEnLossSigmaPrP1);
      fSigDedx[MpdPidUtils::kProton].push_back(fEnLossSigmaPrP2);
      fSigDedx[MpdPidUtils::kProton].push_back(fEnLossSigmaPrP3);
      fdEdxSigmaMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fSigDedx[MpdPidUtils::kProton]));

      TF1 *fEnLossSigmaDeP1 = new TF1("fEnLossSigmaDeP1", "pol1(0)", 0.6, 1.025);
      fEnLossSigmaDeP1->SetParameters(0.158871, -0.0961071);
      TF1 *fEnLossSigmaDeP2 = new TF1("fEnLossSigmaDeP2", "pol1(0)", 1.025, 3.0);
      fEnLossSigmaDeP2->SetParameters(0.054097, 0.00610817);
      TF1 *fEnLossSigmaDeP3 = new TF1("fEnLossSigmaDeP3", "[0]", 3.0, PMAX);
      fEnLossSigmaDeP3->SetParameter(0, fEnLossSigmaDeP2->Eval(3.0));
      fSigDedx[MpdPidUtils::kDeuteron].push_back(fEnLossSigmaDeP1);
      fSigDedx[MpdPidUtils::kDeuteron].push_back(fEnLossSigmaDeP2);
      fSigDedx[MpdPidUtils::kDeuteron].push_back(fEnLossSigmaDeP3);
      fdEdxSigmaMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kDeuteron, fSigDedx[MpdPidUtils::kDeuteron]));

      TF1 *fEnLossSigmaTrP1 = new TF1("fEnLossSigmaTrP1", "pol1(0)", PMIN, 0.625);
      fEnLossSigmaTrP1->SetParameters(0.50625, -0.65);
      TF1 *fEnLossSigmaTrP2 = new TF1("fEnLossSigmaTrP2", "pol1(0)", 0.625, 1.5);
      fEnLossSigmaTrP2->SetParameters(0.137969, -0.05625);
      TF1 *fEnLossSigmaTrP3 = new TF1("fEnLossSigmaTrP3", "pol1(0)", 1.5, 3.0);
      fEnLossSigmaTrP3->SetParameters(0.0510345, 0.0039953);
      TF1 *fEnLossSigmaTrP4 = new TF1("fEnLossSigmaTrP4", "[0]", 3.0, PMAX);
      fEnLossSigmaTrP4->SetParameter(0, fEnLossSigmaTrP3->Eval(3.0));
      fSigDedx[MpdPidUtils::kTriton].push_back(fEnLossSigmaTrP1);
      fSigDedx[MpdPidUtils::kTriton].push_back(fEnLossSigmaTrP2);
      fSigDedx[MpdPidUtils::kTriton].push_back(fEnLossSigmaTrP3);
      fSigDedx[MpdPidUtils::kTriton].push_back(fEnLossSigmaTrP4);
      fdEdxSigmaMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kTriton, fSigDedx[MpdPidUtils::kTriton]));

      TF1 *fEnLossSigmaHe3P1 = new TF1("fEnLossSigmaHe3P1", "pol1(0)", PMIN, 0.675);
      fEnLossSigmaHe3P1->SetParameters(0.44225, -0.47);
      TF1 *fEnLossSigmaHe3P2 = new TF1("fEnLossSigmaHe3P2", "pol1(0)", 0.675, 1.475);
      fEnLossSigmaHe3P2->SetParameters(0.17984375, -0.08125);
      TF1 *fEnLossSigmaHe3P3 = new TF1("fEnLossSigmaHe3P3", "pol1(0)", 1.475, 1.6);
      fEnLossSigmaHe3P3->SetParameters(-0.261818, 0.21818);
      TF1 *fEnLossSigmaHe3P4 = new TF1("fEnLossSigmaHe3P4", "[0]", 1.6, PMAX);
      fEnLossSigmaHe3P4->SetParameter(0, fEnLossSigmaHe3P3->Eval(1.6));
      fSigDedx[MpdPidUtils::kHe3].push_back(fEnLossSigmaHe3P1);
      fSigDedx[MpdPidUtils::kHe3].push_back(fEnLossSigmaHe3P2);
      fSigDedx[MpdPidUtils::kHe3].push_back(fEnLossSigmaHe3P3);
      fSigDedx[MpdPidUtils::kHe3].push_back(fEnLossSigmaHe3P4);
      fdEdxSigmaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe3, fSigDedx[MpdPidUtils::kHe3]));

      TF1 *fEnLossSigmaHe4P1 = new TF1("fEnLossSigmaHe4P1", "pol1(0)", PMIN, 0.475);
      fEnLossSigmaHe4P1->SetParameters(0.37188, -0.225);
      TF1 *fEnLossSigmaHe4P2 = new TF1("fEnLossSigmaHe4P2", "pol1(0)", 0.475, 0.7);
      fEnLossSigmaHe4P2->SetParameters(0.19375, 0.15);
      TF1 *fEnLossSigmaHe4P3 = new TF1("fEnLossSigmaHe4P3", "pol1(0)", 0.7, 0.925);
      fEnLossSigmaHe4P3->SetParameters(0.72231, -0.6025);
      TF1 *fEnLossSigmaHe4P4 = new TF1("fEnLossSigmaHe4P4", "pol1(0)", 0.925, 1.6);
      fEnLossSigmaHe4P4->SetParameters(0.32477, -0.17273);
      TF1 *fEnLossSigmaHe4P5 = new TF1("fEnLossSigmaHe4P5", "[0]", 1.6, PMAX);
      fEnLossSigmaHe4P5->SetParameter(0, fEnLossSigmaHe4P4->Eval(1.6));
      fSigDedx[MpdPidUtils::kHe4].push_back(fEnLossSigmaHe4P1);
      fSigDedx[MpdPidUtils::kHe4].push_back(fEnLossSigmaHe4P2);
      fSigDedx[MpdPidUtils::kHe4].push_back(fEnLossSigmaHe4P3);
      fSigDedx[MpdPidUtils::kHe4].push_back(fEnLossSigmaHe4P4);
      fSigDedx[MpdPidUtils::kHe4].push_back(fEnLossSigmaHe4P5);
      fdEdxSigmaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe4, fSigDedx[MpdPidUtils::kHe4]));

      /// Asymmetry parameter of dE/dx
      TF1 *fAsymElP1 = new TF1("fAsymElP1", "pol1(0)", PMIN, 0.5);
      fAsymElP1->SetParameters(-0.340074, 1.04407);
      TF1 *fAsymElP2 = new TF1("fAsymElP2", "pol1(0)", 0.5, 0.85);
      fAsymElP2->SetParameters(0.291559, -0.306138);
      TF1 *fAsymElP3 = new TF1("fAsymElP3", "pol1(0)", 0.85, 3.0);
      fAsymElP3->SetParameters(0.0915437, 0.0467649);
      TF1 *fAsymElP4 = new TF1("fAsymElP4", "[0]", 3.0, PMAX);
      fAsymElP4->SetParameter(0, fAsymElP3->Eval(3.0));
      fAsymDedx[MpdPidUtils::kElectron].push_back(fAsymElP1);
      fAsymDedx[MpdPidUtils::kElectron].push_back(fAsymElP2);
      fAsymDedx[MpdPidUtils::kElectron].push_back(fAsymElP3);
      fAsymDedx[MpdPidUtils::kElectron].push_back(fAsymElP4);
      fdEdxDeltaMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kElectron, fAsymDedx[MpdPidUtils::kElectron]));

      TF1 *fAsymMuP1 = new TF1("fAsymMuP1", "pol1(0)", PMIN, 0.125);
      fAsymMuP1->SetParameters(2.5656, -14.391);
      TF1 *fAsymMuP2 = new TF1("fAsymMuP2", "pol1(0)", 0.125, 0.5);
      fAsymMuP2->SetParameters(-0.0467746, 0.496758);
      TF1 *fAsymMuP3 = new TF1("fAsymMuP3", "pol1(0)", 0.5, 3.0);
      fAsymMuP3->SetParameters(0.196422, 0.0117276);
      TF1 *fAsymMuP4 = new TF1("fAsymMuP4", "[0]", 3.0, PMAX);
      fAsymMuP4->SetParameter(0, fAsymMuP3->Eval(3.0));
      fAsymDedx[MpdPidUtils::kMuon].push_back(fAsymMuP1);
      fAsymDedx[MpdPidUtils::kMuon].push_back(fAsymMuP2);
      fAsymDedx[MpdPidUtils::kMuon].push_back(fAsymMuP3);
      fAsymDedx[MpdPidUtils::kMuon].push_back(fAsymMuP4);
      fdEdxDeltaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kMuon, fAsymDedx[MpdPidUtils::kMuon]));

      TF1 *fAsymPiP1 = new TF1("fAsymPiP1", "pol1(0)", PMIN, 0.1841);
      fAsymPiP1->SetParameters(4.51001, -23.116);
      TF1 *fAsymPiP2 = new TF1("fAsymPiP2", "pol1(0)", 0.1841, 0.6);
      fAsymPiP2->SetParameters(0.206221, 0.103783);
      TF1 *fAsymPiP3 = new TF1("fAsymPiP3", "pol1(0)", 0.6, 3.0);
      fAsymPiP3->SetParameters(0.204105, 0.0544095);
      TF1 *fAsymPiP4 = new TF1("fAsymPiP4", "[0]", 3.0, PMAX);
      fAsymPiP4->SetParameter(0, fAsymPiP3->Eval(3.0));
      fAsymDedx[MpdPidUtils::kPion].push_back(fAsymPiP1);
      fAsymDedx[MpdPidUtils::kPion].push_back(fAsymPiP2);
      fAsymDedx[MpdPidUtils::kPion].push_back(fAsymPiP3);
      fAsymDedx[MpdPidUtils::kPion].push_back(fAsymPiP4);
      fdEdxDeltaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fAsymDedx[MpdPidUtils::kPion]));

      TF1 *fAsymKaP1 = new TF1("fAsymKaP1", "pol1(0)", PMIN, 0.5);
      fAsymKaP1->SetParameters(2.3588, -4.87114);
      TF1 *fAsymKaP2 = new TF1("fAsymKaP2", "pol1(0)", 0.5, 1.05);
      fAsymKaP2->SetParameters(0.00526557, 0.239378);
      TF1 *fAsymKaP3 = new TF1("fAsymKaP3", "pol1(0)", 1.05, 3.0);
      fAsymKaP3->SetParameters(0.0892569, 0.119403);
      TF1 *fAsymKaP4 = new TF1("fAsymKaP4", "[0]", 3.0, PMAX);
      fAsymKaP4->SetParameter(0, fAsymKaP3->Eval(3.0));
      fAsymDedx[MpdPidUtils::kKaon].push_back(fAsymKaP1);
      fAsymDedx[MpdPidUtils::kKaon].push_back(fAsymKaP2);
      fAsymDedx[MpdPidUtils::kKaon].push_back(fAsymKaP3);
      fAsymDedx[MpdPidUtils::kKaon].push_back(fAsymKaP4);
      fdEdxDeltaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fAsymDedx[MpdPidUtils::kKaon]));

      TF1 *fAsymPrP1 = new TF1("fAsymPrP1", "pol1(0)", PMIN, 0.25);
      fAsymPrP1->SetParameters(1.68519, -2.75816);
      TF1 *fAsymPrP2 = new TF1("fAsymPrP2", "pol1(0)", 0.25, 0.5);
      fAsymPrP2->SetParameters(1.68519, -2.75816);
      TF1 *fAsymPrP3 = new TF1("fAsymPrP3", "pol1(0)", 0.5, 3.0);
      fAsymPrP3->SetParameters(-0.0761382, 0.186669);
      TF1 *fAsymPrP4 = new TF1("fAsymPrP4", "[0]", 3.0, PMAX);
      fAsymPrP4->SetParameter(0, fAsymPrP3->Eval(3.0));
      fAsymDedx[MpdPidUtils::kProton].push_back(fAsymPrP1);
      fAsymDedx[MpdPidUtils::kProton].push_back(fAsymPrP2);
      fAsymDedx[MpdPidUtils::kProton].push_back(fAsymPrP3);
      fAsymDedx[MpdPidUtils::kProton].push_back(fAsymPrP4);
      fdEdxDeltaMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fAsymDedx[MpdPidUtils::kProton]));
   } else {
      /// dE/dx Bethe-Bloch
      it = fdEdxBBMap.find(MpdPidUtils::kElectron);
      it->second[0]->SetParameters(fCoef * (-0.019326681), fCoef * (1.9201386));
      it = fdEdxBBMap.find(MpdPidUtils::kMuon);
      it->second[0]->SetParameters(fCoef * (0.18313862), 5.8679075, -0.0056352172, 1.9553629, 0.90609069);
      it = fdEdxBBMap.find(MpdPidUtils::kPion);
      it->second[0]->SetParameters(fCoef * (6.6409779), -1.1578231, -0.89821579, 0.39576895, 0.0030161715);
      it = fdEdxBBMap.find(MpdPidUtils::kKaon);
      it->second[0]->SetParameters(fCoef * (0.333125), 1.06618, -0.925413, 1.89814, 0.0278704, fCoef * (0.333125),
                                   1.06618, -0.925413, 1.89814, 0.0278704);
      it = fdEdxBBMap.find(MpdPidUtils::kProton);
      it->second[0]->SetParameters(fCoef * (7.2974), 3.56915, 17.0147, 1.01921, -7.75154, fCoef * (0.208992), 1.35048,
                                   -0.980985, 2.01922, 0.0106046);
      // it->second[1]->SetParameter( 0, fCoef*(1.81344) );
      it = fdEdxBBMap.find(MpdPidUtils::kDeuteron);
      it->second[0]->SetParameters(fCoef * (21.479608), 3.3186673, 16.269037, 0.890677, -4.2241652,
                                   fCoef * (0.20764449), 0.48498402, -0.99221832, 2.0873096, 0.0063239652);
      it = fdEdxBBMap.find(MpdPidUtils::kTriton);
      it->second[0]->SetParameters(fCoef * (2.9608171), 7.5586606, 181.76905, 0.94955581, -15.689995,
                                   fCoef * (0.43139727), -1.1708504, -0.97933102, 2.2046764, 0.017628675);
      it = fdEdxBBMap.find(MpdPidUtils::kHe3);
      it->second[0]->SetParameters(fCoef * (-20.795432), -0.47384471, 1.7198872, 1.1896359, fCoef * (-6.6117106),
                                   -0.46247527, 0.33893895, 2.632457);
      it->second[1]->SetParameter(0, it->second[0]->Eval(2.6));
      it = fdEdxBBMap.find(MpdPidUtils::kHe4);
      it->second[0]->SetParameters(fCoef * (-21.022252), -0.49777446, 1.868153, 1.1551161, fCoef * (-2.9013868),
                                   -1.3956671, 0.32487842, 2.1497518);
      it->second[1]->SetParameter(0, it->second[0]->Eval(2.6));
      /// m2
      it = fParM2Map.find(MpdPidUtils::kElectron);
      it->second[0]->SetParameters(0.001227, -0.000973509, 0.0314155);
      it = fParM2Map.find(MpdPidUtils::kMuon);
      it->second[0]->SetParameters(0.00166279, -0.00131341, 0.0311028);
      it = fParM2Map.find(MpdPidUtils::kPion);
      it->second[1]->SetParameters(-0.000480833, 0.00274342, 0.0205474);
      it->second[2]->SetParameters(-0.0307661, 0.0404203, 0.00870571);
      it->second[0]->SetParameter(0, it->second[1]->Eval(it->second[1]->GetXmin()));
      it = fParM2Map.find(MpdPidUtils::kKaon);
      it->second[0]->SetParameters(0.00780637, -0.00395721, 0.029572, -0.00282261);
      it = fParM2Map.find(MpdPidUtils::kProton);
      it->second[1]->SetParameters(0.086827, -0.116522, 0.0874504);
      it->second[2]->SetParameters(0.0240519, 0.0166668, 0.0158857);
      it->second[0]->SetParameter(0, it->second[1]->Eval(it->second[1]->GetXmin()));
      it = fParM2Map.find(MpdPidUtils::kDeuteron);
      it->second[1]->SetParameters(0.200874, -0.163555, 0.0883628);
      it->second[0]->SetParameter(0, it->second[1]->Eval(it->second[1]->GetXmin()));
      it = fParM2Map.find(MpdPidUtils::kTriton);
      it->second[1]->SetParameters(0.948708, -0.474263, 0.151498);
      it->second[0]->SetParameter(0, it->second[1]->Eval(it->second[1]->GetXmin()));
      it = fParM2Map.find(MpdPidUtils::kHe3);
      it->second[1]->SetParameters(0.176375, -0.0917918, 0.0429879);
      it->second[0]->SetParameter(0, it->second[1]->Eval(it->second[1]->GetXmin()));
      it = fParM2Map.find(MpdPidUtils::kHe4);
      it->second[1]->SetParameters(0.379548, -0.195864, 0.061286);
      it->second[0]->SetParameter(0, it->second[1]->Eval(it->second[1]->GetXmin()));

      /// dE/dx width
      TF1 *fEnLossSigmaElP1 = new TF1("fEnLossSigmaElP1", "pol1(0)", PMIN, 1.2667232);
      fEnLossSigmaElP1->SetParameters(0.065546253, -0.00653691);
      TF1 *fEnLossSigmaElP2 = new TF1("fEnLossSigmaElP2", "pol1(0)", 1.2667232, PMAX);
      fEnLossSigmaElP2->SetParameters(0.056658487, 0.00047943452);
      fSigDedx[MpdPidUtils::kElectron].push_back(fEnLossSigmaElP1);
      fSigDedx[MpdPidUtils::kElectron].push_back(fEnLossSigmaElP2);
      fdEdxSigmaMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kElectron, fSigDedx[MpdPidUtils::kElectron]));

      TF1 *fEnLossSigmaMuP1 = new TF1("fEnLossSigmaMuP1", "pol1(0)", PMIN, 0.23238453);
      fEnLossSigmaMuP1->SetParameters(0.096834625, -0.088311039);
      TF1 *fEnLossSigmaMuP2 = new TF1("fEnLossSigmaMuP2", "pol1(0)", 0.23238453, 1.0873775);
      fEnLossSigmaMuP2->SetParameters(0.078958282, -0.011385339);
      TF1 *fEnLossSigmaMuP3 = new TF1("fEnLossSigmaMuP3", "pol1(0)", 1.0873775, 1.7854915);
      fEnLossSigmaMuP3->SetParameters(0.070961629, -0.004031266);
      TF1 *fEnLossSigmaMuP4 = new TF1("fEnLossSigmaMuP4", "pol1(0)", 1.7854915, PMAX);
      fEnLossSigmaMuP4->SetParameters(0.065691429, -0.0010795862);
      fSigDedx[MpdPidUtils::kMuon].push_back(fEnLossSigmaMuP1);
      fSigDedx[MpdPidUtils::kMuon].push_back(fEnLossSigmaMuP2);
      fSigDedx[MpdPidUtils::kMuon].push_back(fEnLossSigmaMuP3);
      fSigDedx[MpdPidUtils::kMuon].push_back(fEnLossSigmaMuP4);
      fdEdxSigmaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kMuon, fSigDedx[MpdPidUtils::kMuon]));

      TF1 *fEnLossSigmaPiP1 = new TF1("fEnLossSigmaPiP1", "pol1(0)", PMIN, 0.135226);
      fEnLossSigmaPiP1->SetParameters(0.0978523, -0.0468431);
      TF1 *fEnLossSigmaPiP2 = new TF1("fEnLossSigmaPiP2", "pol1(0)", 0.135226, 0.270254);
      fEnLossSigmaPiP2->SetParameters(0.112118, -0.152342);
      TF1 *fEnLossSigmaPiP3 = new TF1("fEnLossSigmaPiP3", "pol1(0)", 0.270254, 1.96711);
      fEnLossSigmaPiP3->SetParameters(0.0717427, -0.00294257);
      TF1 *fEnLossSigmaPiP4 = new TF1("fEnLossSigmaPiP4", "pol1(0)", 1.96711, 3.0);
      fEnLossSigmaPiP4->SetParameters(0.0567998, 0.00465381);
      TF1 *fEnLossSigmaPiP5 = new TF1("fEnLossSigmaPiP5", "[0]", 3.0, PMAX);
      fEnLossSigmaPiP5->SetParameter(0, fEnLossSigmaPiP4->Eval(3.0));
      fSigDedx[MpdPidUtils::kPion].push_back(fEnLossSigmaPiP1);
      fSigDedx[MpdPidUtils::kPion].push_back(fEnLossSigmaPiP2);
      fSigDedx[MpdPidUtils::kPion].push_back(fEnLossSigmaPiP3);
      fSigDedx[MpdPidUtils::kPion].push_back(fEnLossSigmaPiP4);
      fSigDedx[MpdPidUtils::kPion].push_back(fEnLossSigmaPiP5);
      fdEdxSigmaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fSigDedx[MpdPidUtils::kPion]));

      TF1 *fEnLossSigmaKaP1 = new TF1("fEnLossSigmaKaP1", "pol1(0)", PMIN, 0.52439);
      fEnLossSigmaKaP1->SetParameters(0.114346, -0.0814144);
      TF1 *fEnLossSigmaKaP2 = new TF1("fEnLossSigmaKaP2", "pol1(0)", 0.52439, 3.0);
      fEnLossSigmaKaP2->SetParameters(0.0728159, -0.00221717);
      TF1 *fEnLossSigmaKaP3 = new TF1("fEnLossSigmaKaP3", "[0]", 3.0, PMAX);
      fEnLossSigmaKaP3->SetParameter(0, fEnLossSigmaKaP2->Eval(3.0));
      fSigDedx[MpdPidUtils::kKaon].push_back(fEnLossSigmaKaP1);
      fSigDedx[MpdPidUtils::kKaon].push_back(fEnLossSigmaKaP2);
      fSigDedx[MpdPidUtils::kKaon].push_back(fEnLossSigmaKaP3);
      fdEdxSigmaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fSigDedx[MpdPidUtils::kKaon]));

      TF1 *fEnLossSigmaPrP1 = new TF1("fEnLossSigmaPrP1", "[0]", PMIN, 0.1);
      fEnLossSigmaPrP1->SetParameters(0, 0.21297);
      TF1 *fEnLossSigmaPrP2 = new TF1("fEnLossSigmaPrP2", "pol1(0)", 0.1, 0.246866);
      fEnLossSigmaPrP2->SetParameters(0.291244, -0.782737);
      TF1 *fEnLossSigmaPrP3 = new TF1("fEnLossSigmaPrP3", "pol1(0)", 0.246866, 0.903752);
      fEnLossSigmaPrP3->SetParameters(0.108798, -0.0436887);
      TF1 *fEnLossSigmaPrP4 = new TF1("fEnLossSigmaPrP4", "pol1(0)", 0.903752, 3.0);
      fEnLossSigmaPrP4->SetParameters(0.0714059, -0.00231462);
      TF1 *fEnLossSigmaPrP5 = new TF1("fEnLossSigmaPrP5", "[0]", 3.0, PMAX);
      fEnLossSigmaPrP5->SetParameter(0, fEnLossSigmaPrP4->Eval(3.0));
      fSigDedx[MpdPidUtils::kProton].push_back(fEnLossSigmaPrP1);
      fSigDedx[MpdPidUtils::kProton].push_back(fEnLossSigmaPrP2);
      fSigDedx[MpdPidUtils::kProton].push_back(fEnLossSigmaPrP3);
      fSigDedx[MpdPidUtils::kProton].push_back(fEnLossSigmaPrP4);
      fSigDedx[MpdPidUtils::kProton].push_back(fEnLossSigmaPrP5);
      fdEdxSigmaMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fSigDedx[MpdPidUtils::kProton]));

      TF1 *fEnLossSigmaDeP2 = new TF1("fEnLossSigmaDeP2", "pol1(0)", 0.76311362, 1.1312292);
      fEnLossSigmaDeP2->SetParameters(0.091285714, -0.017142857);
      TF1 *fEnLossSigmaDeP1 = new TF1("fEnLossSigmaDeP1", "[0]", PMIN, 0.76311362);
      fEnLossSigmaDeP1->SetParameter(0, fEnLossSigmaDeP2->Eval(fEnLossSigmaDeP2->GetXmin()));
      TF1 *fEnLossSigmaDeP3 = new TF1("fEnLossSigmaDeP3", "pol1(0)", 1.1312292, 3.0);
      fEnLossSigmaDeP3->SetParameters(0.070917901, 0.00086217052);
      TF1 *fEnLossSigmaDeP4 = new TF1("fEnLossSigmaDeP4", "[0]", 3.0, PMAX);
      fEnLossSigmaDeP4->SetParameter(0, fEnLossSigmaDeP3->Eval(3.0));
      fSigDedx[MpdPidUtils::kDeuteron].push_back(fEnLossSigmaDeP1);
      fSigDedx[MpdPidUtils::kDeuteron].push_back(fEnLossSigmaDeP2);
      fSigDedx[MpdPidUtils::kDeuteron].push_back(fEnLossSigmaDeP3);
      fSigDedx[MpdPidUtils::kDeuteron].push_back(fEnLossSigmaDeP4);
      fdEdxSigmaMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kDeuteron, fSigDedx[MpdPidUtils::kDeuteron]));

      TF1 *fEnLossSigmaTrP1 = new TF1("fEnLossSigmaTrP1", "pol1(0)", PMIN, 3.0);
      fEnLossSigmaTrP1->SetParameters(0.079593779, -0.0016492155);
      TF1 *fEnLossSigmaTrP2 = new TF1("fEnLossSigmaTrP2", "[0]", 3.0, PMAX);
      fEnLossSigmaTrP2->SetParameter(0, fEnLossSigmaTrP1->Eval(3.0));
      fSigDedx[MpdPidUtils::kTriton].push_back(fEnLossSigmaTrP1);
      fSigDedx[MpdPidUtils::kTriton].push_back(fEnLossSigmaTrP2);
      fdEdxSigmaMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kTriton, fSigDedx[MpdPidUtils::kTriton]));

      TF1 *fEnLossSigmaHe3P1 = new TF1("fEnLossSigmaHe3P1", "pol1(0)", PMIN, 1.5222468);
      fEnLossSigmaHe3P1->SetParameters(0.084763216, -0.022888375);
      TF1 *fEnLossSigmaHe3P2 = new TF1("fEnLossSigmaHe3P2", "pol1(0)", 1.5222468, 2.6);
      fEnLossSigmaHe3P2->SetParameters(0.059599353, -0.006357636);
      TF1 *fEnLossSigmaHe3P3 = new TF1("fEnLossSigmaHe3P3", "[0]", 2.6, PMAX);
      fEnLossSigmaHe3P3->SetParameter(0, fEnLossSigmaHe3P2->Eval(2.6));
      fSigDedx[MpdPidUtils::kHe3].push_back(fEnLossSigmaHe3P1);
      fSigDedx[MpdPidUtils::kHe3].push_back(fEnLossSigmaHe3P2);
      fSigDedx[MpdPidUtils::kHe3].push_back(fEnLossSigmaHe3P3);
      fdEdxSigmaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe3, fSigDedx[MpdPidUtils::kHe3]));

      TF1 *fEnLossSigmaHe4P1 = new TF1("fEnLossSigmaHe4P1", "pol1(0)", PMIN, 2.6);
      fEnLossSigmaHe4P1->SetParameters(0.067640377, -0.0090089899);
      TF1 *fEnLossSigmaHe4P2 = new TF1("fEnLossSigmaHe4P2", "[0]", 2.6, PMAX);
      fEnLossSigmaHe4P2->SetParameter(0, fEnLossSigmaHe4P1->Eval(2.6));
      fSigDedx[MpdPidUtils::kHe4].push_back(fEnLossSigmaHe4P1);
      fSigDedx[MpdPidUtils::kHe4].push_back(fEnLossSigmaHe4P2);
      fdEdxSigmaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe4, fSigDedx[MpdPidUtils::kHe4]));

      /// Asymmetry parameter of dE/dx
      TF1 *fAsymElP1 = new TF1("fAsymElP1", "[0]", PMIN, 2.3676103);
      fAsymElP1->SetParameter(0, 0.5);
      TF1 *fAsymElP2 = new TF1("fAsymElP2", "pol1(0)", 2.3676103, PMAX);
      fAsymElP2->SetParameters(0.68870097, -0.079701026);
      fAsymDedx[MpdPidUtils::kElectron].push_back(fAsymElP1);
      fAsymDedx[MpdPidUtils::kElectron].push_back(fAsymElP2);
      fdEdxDeltaMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kElectron, fAsymDedx[MpdPidUtils::kElectron]));

      TF1 *fAsymMuP1 = new TF1("fAsymMuP1", "pol1(0)", PMIN, 0.17693005);
      fAsymMuP1->SetParameters(4.3617776, -21.934439);
      TF1 *fAsymMuP2 = new TF1("fAsymMuP2", "pol1(0)", 0.17693005, 1.0495326);
      fAsymMuP2->SetParameters(0.47825147, 0.015061837);
      TF1 *fAsymMuP3 = new TF1("fAsymMuP3", "pol1(0)", 1.0495326, PMAX);
      fAsymMuP3->SetParameters(0.54616505, -0.049646561);
      fAsymDedx[MpdPidUtils::kMuon].push_back(fAsymMuP1);
      fAsymDedx[MpdPidUtils::kMuon].push_back(fAsymMuP2);
      fAsymDedx[MpdPidUtils::kMuon].push_back(fAsymMuP3);
      fdEdxDeltaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kMuon, fAsymDedx[MpdPidUtils::kMuon]));

      TF1 *fAsymPiP1 = new TF1("fAsymPiP1", "pol1(0)", PMIN, 0.183903);
      fAsymPiP1->SetParameters(4.573, -22.8192);
      TF1 *fAsymPiP2 = new TF1("fAsymPiP2", "pol1(0)", 0.183903, 0.559408);
      fAsymPiP2->SetParameters(0.405159, -0.155943);
      TF1 *fAsymPiP3 = new TF1("fAsymPiP3", "pol1(0)", 0.559408, 3.0);
      fAsymPiP3->SetParameters(0.351838, -0.0606272);
      TF1 *fAsymPiP4 = new TF1("fAsymPiP4", "[0]", 3.0, PMAX);
      fAsymPiP4->SetParameter(0, fAsymPiP3->Eval(3.0));
      fAsymDedx[MpdPidUtils::kPion].push_back(fAsymPiP1);
      fAsymDedx[MpdPidUtils::kPion].push_back(fAsymPiP2);
      fAsymDedx[MpdPidUtils::kPion].push_back(fAsymPiP3);
      fAsymDedx[MpdPidUtils::kPion].push_back(fAsymPiP4);
      fdEdxDeltaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fAsymDedx[MpdPidUtils::kPion]));

      TF1 *fAsymKaP1 = new TF1("fAsymKaP1", "[0]", PMIN, 0.176873);
      fAsymKaP1->SetParameters(0, 2.0);
      TF1 *fAsymKaP2 = new TF1("fAsymKaP2", "pol1(0)", 0.176873, 0.292257);
      fAsymKaP2->SetParameters(4.42093, -13.6874);
      TF1 *fAsymKaP3 = new TF1("fAsymKaP3", "pol1(0)", 0.292257, 0.457227);
      fAsymKaP3->SetParameters(0.623815, -0.694969);
      TF1 *fAsymKaP4 = new TF1("fAsymKaP4", "pol1(0)", 0.457227, 3.0);
      fAsymKaP4->SetParameters(0.296643, 0.0205892);
      TF1 *fAsymKaP5 = new TF1("fAsymKaP5", "[0]", 3.0, PMAX);
      fAsymKaP5->SetParameter(0, fAsymKaP4->Eval(3.0));
      fAsymDedx[MpdPidUtils::kKaon].push_back(fAsymKaP1);
      fAsymDedx[MpdPidUtils::kKaon].push_back(fAsymKaP2);
      fAsymDedx[MpdPidUtils::kKaon].push_back(fAsymKaP3);
      fAsymDedx[MpdPidUtils::kKaon].push_back(fAsymKaP4);
      fAsymDedx[MpdPidUtils::kKaon].push_back(fAsymKaP5);
      fdEdxDeltaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fAsymDedx[MpdPidUtils::kKaon]));

      TF1 *fAsymPrP1 = new TF1("fAsymPrP1", "[0]", PMIN, 0.2);
      fAsymPrP1->SetParameter(0, 0.858727); /// <...>P2.Eval(0.15958801)
      TF1 *fAsymPrP2 = new TF1("fAsymPrP2", "pol1(0)", 0.2, 0.487588);
      fAsymPrP2->SetParameters(1.3097, -2.25487);
      TF1 *fAsymPrP3 = new TF1("fAsymPrP3", "pol1(0)", 0.487588, 0.664107);
      fAsymPrP3->SetParameters(0.264705, -0.111677);
      TF1 *fAsymPrP4 = new TF1("fAsymPrP4", "pol1(0)", 0.664107, 1.38463);
      fAsymPrP4->SetParameters(0.0612819, 0.194634);
      TF1 *fAsymPrP5 = new TF1("fAsymPrP5", "pol1(0)", 1.38463, 3.0);
      fAsymPrP5->SetParameters(0.253465, 0.0558374);
      TF1 *fAsymPrP6 = new TF1("fAsymPrP6", "[0]", 3.0, PMAX);
      fAsymPrP6->SetParameter(0, fAsymPrP4->Eval(3.0));
      fAsymDedx[MpdPidUtils::kProton].push_back(fAsymPrP1);
      fAsymDedx[MpdPidUtils::kProton].push_back(fAsymPrP2);
      fAsymDedx[MpdPidUtils::kProton].push_back(fAsymPrP3);
      fAsymDedx[MpdPidUtils::kProton].push_back(fAsymPrP4);
      fAsymDedx[MpdPidUtils::kProton].push_back(fAsymPrP5);
      fAsymDedx[MpdPidUtils::kProton].push_back(fAsymPrP6);
      fdEdxDeltaMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fAsymDedx[MpdPidUtils::kProton]));

      TF1 *fAsymDeP1 = new TF1("fAsymDeP1", "pol1(0)", PMIN, 0.46066487);
      fAsymDeP1->SetParameters(-0.12857143, 1.3714286);
      TF1 *fAsymDeP2 = new TF1("fAsymDeP2", "pol1(0)", 0.46066487, 3.0);
      fAsymDeP2->SetParameters(0.4816788, 0.046712329);
      TF1 *fAsymDeP3 = new TF1("fAsymDeP3", "[0]", 3.0, PMAX);
      fAsymDeP3->SetParameter(0, fAsymDeP2->Eval(3.0));
      fAsymDedx[MpdPidUtils::kDeuteron].push_back(fAsymDeP1);
      fAsymDedx[MpdPidUtils::kDeuteron].push_back(fAsymDeP2);
      fAsymDedx[MpdPidUtils::kDeuteron].push_back(fAsymDeP3);
      fdEdxDeltaMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kDeuteron, fAsymDedx[MpdPidUtils::kDeuteron]));

      TF1 *fAsymTrP1 = new TF1("fAsymTrP1", "pol1(0)", PMIN, 0.83782809);
      fAsymTrP1->SetParameters(-0.1027975, 0.79037208);
      TF1 *fAsymTrP2 = new TF1("fAsymTrP2", "[0]", 0.83782809, PMAX);
      fAsymTrP2->SetParameter(0, 0.55939843);
      fAsymDedx[MpdPidUtils::kTriton].push_back(fAsymTrP1);
      fAsymDedx[MpdPidUtils::kTriton].push_back(fAsymTrP2);
      fdEdxDeltaMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kTriton, fAsymDedx[MpdPidUtils::kTriton]));

      TF1 *fAsymHe3P1 = new TF1("fAsymHe3P1", "pol1(0)", PMIN, 0.47973371);
      fAsymHe3P1->SetParameters(-0.51177512, 0.7596956);
      TF1 *fAsymHe3P2 = new TF1("fAsymHe3P2", "pol1(0)", 0.47973371, 1.0043019);
      fAsymHe3P2->SetParameters(-0.67019153, 1.089913);
      TF1 *fAsymHe3P3 = new TF1("fAsymHe3P3", "pol1(0)", 1.0043019, 2.6);
      fAsymHe3P3->SetParameters(0.33736312, 0.086674133);
      TF1 *fAsymHe3P4 = new TF1("fAsymHe3P4", "[0]", 2.6, PMAX);
      fAsymHe3P4->SetParameter(0, fAsymHe3P3->Eval(2.6));
      fAsymDedx[MpdPidUtils::kHe3].push_back(fAsymHe3P1);
      fAsymDedx[MpdPidUtils::kHe3].push_back(fAsymHe3P2);
      fAsymDedx[MpdPidUtils::kHe3].push_back(fAsymHe3P3);
      fAsymDedx[MpdPidUtils::kHe3].push_back(fAsymHe3P4);
      fdEdxDeltaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe3, fAsymDedx[MpdPidUtils::kHe3]));

      TF1 *fAsymHe4P1 = new TF1("fAsymHe4P1", "pol1(0)", PMIN, 0.6096774);
      fAsymHe4P1->SetParameters(-0.36658037, 0.18506071);
      TF1 *fAsymHe4P2 = new TF1("fAsymHe4P2", "pol1(0)", 0.6096774, 1.0146462);
      fAsymHe4P2->SetParameters(-0.91671437, 1.087397);
      TF1 *fAsymHe4P3 = new TF1("fAsymHe4P3", "pol1(0)", 1.0146462, 1.4051261);
      fAsymHe4P3->SetParameters(-0.68765897, 0.86164797);
      TF1 *fAsymHe4P4 = new TF1("fAsymHe4P4", "pol1(0)", 1.4051261, 2.6);
      fAsymHe4P4->SetParameters(0.53212178, -0.0064454674);
      TF1 *fAsymHe4P5 = new TF1("fAsymHe4P5", "[0]", 2.6, PMAX);
      fAsymHe4P5->SetParameter(0, fAsymHe4P4->Eval(2.6));
      fAsymDedx[MpdPidUtils::kHe4].push_back(fAsymHe4P1);
      fAsymDedx[MpdPidUtils::kHe4].push_back(fAsymHe4P2);
      fAsymDedx[MpdPidUtils::kHe4].push_back(fAsymHe4P3);
      fAsymDedx[MpdPidUtils::kHe4].push_back(fAsymHe4P4);
      fAsymDedx[MpdPidUtils::kHe4].push_back(fAsymHe4P5);
      fdEdxDeltaMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe4, fAsymDedx[MpdPidUtils::kHe4]));
   }

   /// Particle yields versus momentum (close to a thermal function).
   /// The predicted number roughly speaking are not dN/dy or total yiels,
   /// only their relative hights are of relevance!

   Double_t amplParam;

   TF1 *fParElPosMom = new TF1("fParElPosMom", this, &MpdPid::MomPr, PMIN, PMAX, 5, "MpdPid", "MomPi");
   TF1 *fParElNegMom = new TF1("fParElNegMom", this, &MpdPid::MomPr, PMIN, PMAX, 5, "MpdPid", "MomPi");
   TF1 *fParMuPosMom = new TF1("fParMuPosMom", this, &MpdPid::MomPi, PMIN, PMAX, 5, "MpdPid", "MomPi");
   TF1 *fParMuNegMom = new TF1("fParMuNegMom", this, &MpdPid::MomPi, PMIN, PMAX, 5, "MpdPid", "MomPi");
   TF1 *fParPiPosMom = new TF1("fParPiPosMom", this, &MpdPid::MomPi, PMIN, PMAX, 5, "MpdPid", "MomPi");
   TF1 *fParPiNegMom = new TF1("fParPiNegMom", this, &MpdPid::MomPi, PMIN, PMAX, 5, "MpdPid", "MomPi");
   TF1 *fParKaPosMom = new TF1("fParKaPosMom", this, &MpdPid::MomPr, PMIN, PMAX, 5, "MpdPid", "MomPr");
   TF1 *fParKaNegMom = new TF1("fParKaNegMom", this, &MpdPid::MomPr, PMIN, PMAX, 5, "MpdPid", "MomPr");
   TF1 *fParPrPosMom = new TF1("fParPrPosMom", this, &MpdPid::MomPr, PMIN, PMAX, 5, "MpdPid", "MomPr");
   TF1 *fParPrNegMom = new TF1("fParPrNegMom", this, &MpdPid::MomPr, PMIN, PMAX, 5, "MpdPid", "MomPr");
   TF1 *fParDeMom    = new TF1("fParDeMom", this, &MpdPid::MomPr, PMIN, PMAX, 5, "MpdPid", "MomPr");
   TF1 *fParTrMom    = new TF1("fParTrMom", this, &MpdPid::MomPr, PMIN, PMAX, 5, "MpdPid", "MomPr");
   TF1 *fParHe3Mom   = new TF1("fParHe3Mom", this, &MpdPid::MomPr, PMIN, PMAX, 5, "MpdPid", "MomPr");
   TF1 *fParHe4Mom   = new TF1("fParHe4Mom", this, &MpdPid::MomPr, PMIN, PMAX, 5, "MpdPid", "MomPr");

   if (Generator == "NSIG") {
      fMethod = kFALSE;
      delete fParElPosMom;
      delete fParElNegMom;
      delete fParMuPosMom;
      delete fParMuNegMom;
      delete fParPiPosMom;
      delete fParPiNegMom;
      delete fParKaPosMom;
      delete fParKaNegMom;
      delete fParPrPosMom;
      delete fParPrNegMom;
      delete fParDeMom;
      delete fParTrMom;
      delete fParHe3Mom;
      delete fParHe4Mom;
      if (NSigPart.Contains("el")) {
         fNSigSpecies[MpdPidUtils::kElectron] = kTRUE;
         cout << "electrons are included in n-sigma method..." << endl;
      } else {
         fNSigSpecies[MpdPidUtils::kElectron] = kFALSE;
      }
      if (NSigPart.Contains("mu")) {
         fNSigSpecies[MpdPidUtils::kMuon] = kTRUE;
         cout << "muons are included in n-sigma method..." << endl;
      } else {
         fNSigSpecies[MpdPidUtils::kMuon] = kFALSE;
      }
      if (NSigPart.Contains("pi")) {
         fNSigSpecies[MpdPidUtils::kPion] = kTRUE;
         cout << "pions are included in n-sigma method..." << endl;
      } else {
         fNSigSpecies[MpdPidUtils::kPion] = kFALSE;
      }
      if (NSigPart.Contains("ka")) {
         fNSigSpecies[MpdPidUtils::kKaon] = kTRUE;
         cout << "kaons are included in n-sigma method..." << endl;
      } else {
         fNSigSpecies[MpdPidUtils::kKaon] = kFALSE;
      }
      if (NSigPart.Contains("pr")) {
         fNSigSpecies[MpdPidUtils::kProton] = kTRUE;
         cout << "(anti-)protons are included in n-sigma method..." << endl;
      } else {
         fNSigSpecies[MpdPidUtils::kProton] = kFALSE;
      }
      if (NSigPart.Contains("de")) {
         fNSigSpecies[MpdPidUtils::kDeuteron] = kTRUE;
         cout << "deuterons are included in n-sigma method..." << endl;
      } else {
         fNSigSpecies[MpdPidUtils::kDeuteron] = kFALSE;
      }
      if (NSigPart.Contains("tr")) {
         fNSigSpecies[MpdPidUtils::kTriton] = kTRUE;
         cout << "tritons are included in n-sigma method..." << endl;
      } else {
         fNSigSpecies[MpdPidUtils::kTriton] = kFALSE;
      }
      if (NSigPart.Contains("he3")) {
         fNSigSpecies[MpdPidUtils::kHe3] = kTRUE;
         cout << "he3 are included in n-sigma method..." << endl;
      } else {
         fNSigSpecies[MpdPidUtils::kHe3] = kFALSE;
      }
      if (NSigPart.Contains("he4")) {
         fNSigSpecies[MpdPidUtils::kHe4] = kTRUE;
         cout << "he4 are included in n-sigma method..." << endl;
      } else {
         fNSigSpecies[MpdPidUtils::kHe4] = kFALSE;
      }
   } else
      fMethod = kTRUE;

   if (!((Generator == "LAQGSM") || (Generator == "QGSM") || (Generator == "URQMD") || (Generator == "NSIG") ||
         (Generator == "PHSD") || (Generator == "EPOS") || (Generator == "PHSD_CENT") || (Generator == "PHSD_CSR") ||
         (Generator == "PHSD_NOCSR") || (Generator == "PHQMD"))) {
      cout << "Incorrect generator string! Switching to DEFAULT..." << endl;
      Generator = "DEFAULT";
   }

   if (Generator == "REQ25") {
      fParPiPosMom->SetParameters(315296, 1.26945, 0.155599, 0.331509, 0.238229);
      fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
      fParPiNegMom->SetParameters(45212.8, 7.02777, 0.348156, 0.17387, -0.153182);
      fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
      fPartYieldMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
      fParKaPosMom->SetParameters(111833, 0.27357, 0.347375, 0.11195, 0.30159);
      fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
      fParKaNegMom->SetParameters(62516, 0.252819, 0.301713, 0.0847874, 0.382075);
      fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
      fPartYieldMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
      fParPrPosMom->SetParameters(588372, 1.30255, 0.372303, 0.250823, 0.354083);
      fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
      fParPrNegMom->SetParameters(6604.83, 0.594007, 0.301938, 0.13116, 0.992751);
      fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
      fPartYieldMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
   }

   if (Generator == "EPOS") { /// for p + p collisions
      if (fEnergy < 7.0) {    /// sqrt(s) = 6 GeV, ~1M events
         fParPiPosMom->SetParameters(2853.8, 0.623542, 0.0537214, 0.137323, 1.4267);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
         fParPiNegMom->SetParameters(1407.44, 0.987421, 0.0574761, 0.131526, 1.18167);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
         fParKaPosMom->SetParameters(228.358, 7.59862, 0.162463, 0.148098, 0.183593);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
         fParKaNegMom->SetParameters(57.3598, 0.546853, 0.154656, 0.0593315, 0.543194);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
         fParPrPosMom->SetParameters(1796.08, 8.20092, 0.272697, 0.253002, 0.181719);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
         fParPrNegMom->SetParameters(38.0024, 0.0768751, 0.151354, 0.144691, -0.144209);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
      } else if (fEnergy < 9.0) { /// sqrt(s) = 8.76 GeV, 1M events
         fParPiPosMom->SetParameters(3160.76, 0.70523, 0.080124, 0.206086, 0.974646);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
         fParPiNegMom->SetParameters(1751.43, 1.17603, 0.0914196, 0.206213, 0.737364);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
         fParKaPosMom->SetParameters(362.313, 0.332122, 0.232739, 0.0763912, 0.301778);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
         fParKaNegMom->SetParameters(161.454, 0.483719, 0.217565, 0.0874472, 0.345383);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
         fParPrPosMom->SetParameters(1254.12, 0.409612, 0.385691, 0.171431, 0.181715);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
         fParPrNegMom->SetParameters(89.4667, 0.216063, 0.139082, 0.0149435, 1.29752);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
      } else if (fEnergy < 12.0) { /// sqrt(s) = 10 GeV, 1M events
         fParPiPosMom->SetParameters(3320.92, 0.698313, 0.0892079, 0.227947, 0.875591);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
         fParPiNegMom->SetParameters(2090.76, 1.0203, 0.100926, 0.228817, 0.673615);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
         fParKaPosMom->SetParameters(400.101, 4.6224, 0.208173, 0.182652, 0.069632);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
         fParKaNegMom->SetParameters(197.395, 0.355122, 0.235916, 0.0794332, 0.358637);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
         fParPrPosMom->SetParameters(1089.56, 0.242752, 0.407993, 0.13063, 0.250824);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
         fParPrNegMom->SetParameters(136.99, 1.09482, 0.231846, 0.127863, 0.205208);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
      } else if (fEnergy < 16.0) { /// sqrt(s) = 15 GeV, 1M events
         fParPiPosMom->SetParameters(3750.67, 0.657879, 0.108281, 0.279038, 0.725381);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
         fParPiNegMom->SetParameters(2984.16, 0.769168, 0.129349, 0.287196, 0.507861);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
         fParKaPosMom->SetParameters(675.506, 0.0596585, 0.26442, 0.14352, -0.138467);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
         fParKaNegMom->SetParameters(292.368, 0.256861, 0.275726, 0.0776192, 0.326165);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
         fParPrPosMom->SetParameters(890.967, 0.00141238, 0.416039, 0.150785, -0.150519);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
         fParPrNegMom->SetParameters(518.922, 0.00427706, 0.2737, 0.156507, -0.156194);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
      } else if (fEnergy < 19.0) { /// sqrt(s) = 17.3 GeV, 1M events
         fParPiPosMom->SetParameters(3817.06, 0.664315, 0.115249, 0.292103, 0.674843);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
         fParPiNegMom->SetParameters(3440.51, 0.65109, 0.139433, 0.3066, 0.464094);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
         fParKaPosMom->SetParameters(505.915, 0.229333, 0.289876, 0.0773325, 0.275986);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
         fParKaNegMom->SetParameters(321.28, 0.255435, 0.283842, 0.0806211, 0.315722);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
         fParPrPosMom->SetParameters(655.647, 0.165912, 0.410953, 0.0981928, 0.296843);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
         fParPrNegMom->SetParameters(368.293, 0.275794, 0.289646, 0.066954, 0.501408);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
      } else if (fEnergy < 23.0) { /// sqrt(s) = 20 GeV, 1M events
         fParPiPosMom->SetParameters(3908.67, 0.663881, 0.11984, 0.303021, 0.648025);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
         fParPiNegMom->SetParameters(3652.92, 0.635229, 0.145063, 0.317634, 0.434774);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
         fParKaPosMom->SetParameters(1504.68, -0.608577, 0.300936, 0.105427, -0.169068);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
         fParKaNegMom->SetParameters(352.373, 0.234724, 0.293469, 0.0797531, 0.303124);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
         fParPrPosMom->SetParameters(599.226, 0.172265, 0.406555, 0.0998964, 0.284495);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
         fParPrNegMom->SetParameters(441.899, 0.187852, 0.296744, 0.0391455, 0.591304);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
      } else { /// sqrt(s) = 25 GeV, 1M events
         fParPiPosMom->SetParameters(3982.71, 0.68892, 0.122378, 0.312125, 0.636517);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
         fParPiNegMom->SetParameters(3980.97, 0.60921, 0.149336, 0.331234, 0.424032);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
         fParKaPosMom->SetParameters(1370.93, -0.538074, 0.313878, 0.103715, -0.16116);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
         fParKaNegMom->SetParameters(396.416, 0.206327, 0.304933, 0.0762121, 0.305453);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
         fParPrPosMom->SetParameters(547.27, 0.166229, 0.401412, 0.0955138, 0.284058);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
         fParPrNegMom->SetParameters(546.003, 0.664169, 0.32501, 0.154779, -0.0269533);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
      }
   }
   if (Generator == "PHQMD") { /// Mass Production, fTracking = 2, 15M min bias events @ 8.8 GeV
      fParPiPosMom->SetParameters(245004., 6.43066, 0.285935, 0.191311, 0.0632783);
      fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
      fParPiNegMom->SetParameters(63866.2, 29.1281, 0.284521, 0.203969, 0.0623634);
      fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
      fPartYieldMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
      fParKaPosMom->SetParameters(410468., 0.431534, 0.276289, 0.108674, 0.370743);
      fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
      fParKaNegMom->SetParameters(178284., 0.426688, 0.257688, 0.0968996, 0.439992);
      fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
      fPartYieldMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
      fParPrPosMom->SetParameters(949875., 3.6411, 0.314117, 0.264786, 0.365604);
      fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
      fParPrNegMom->SetParameters(10839.6, 0.705385, 0.221633, 0.0956065, 1.45505);
      fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
      fPartYieldMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
      fParDeMom->SetParameters(15978.8, 2.26702, 0.579366, 0.441819, 0.971303);
      fPartYield[MpdPidUtils::kDeuteron].push_back(fParDeMom);
      fPartYieldMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kDeuteron, fPartYield[MpdPidUtils::kDeuteron]));
      fParTrMom->SetParameters(408.611, 0.184547, 0.589216, 0.103737, 7.26799);
      fPartYield[MpdPidUtils::kTriton].push_back(fParTrMom);
      fPartYieldMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kTriton, fPartYield[MpdPidUtils::kTriton]));
      fParHe3Mom->SetParameters(379.097, 5.99965, 0.327897, 0.289616, 1.33579);
      fPartYield[MpdPidUtils::kHe3].push_back(fParHe3Mom);
      fPartYieldMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe3, fPartYield[MpdPidUtils::kHe3]));
      fParHe4Mom->SetParameters(14.4574, 12.1212, 0.0348463, 0.0324572, 33.2308);
      fPartYield[MpdPidUtils::kHe4].push_back(fParHe4Mom);
      fPartYieldMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe4, fPartYield[MpdPidUtils::kHe4]));
   }
   if (Generator == "PHSD") {
      if (fTrackingState == MpdPidUtils::kHP) {
         if (fEnergy < 7.0) { /// not ready, QGSM 5 gev
            fParElPosMom->SetParameters(17.6, -0.12, 0.078, 0.167, 0.00);
            fPartYield[MpdPidUtils::kElectron].push_back(fParElPosMom);
            fParElNegMom->SetParameters(16.3, -0.12, 0.078, 0.167, 0.00);
            fPartYield[MpdPidUtils::kElectron].push_back(fParElNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kElectron, fPartYield[MpdPidUtils::kElectron]));
            fParMuPosMom->SetParameters(20.5, 0.064, 0.107, 0.05, 0.105);
            fPartYield[MpdPidUtils::kMuon].push_back(fParMuPosMom);
            fParMuNegMom->SetParameters(20.5, 0.064, 0.107, 0.05, 0.105);
            fPartYield[MpdPidUtils::kMuon].push_back(fParMuNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kMuon, fPartYield[MpdPidUtils::kMuon]));
            fParPiPosMom->SetParameters(307.0, 0.035, 0.175, 0.127, 0.139);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
            fParPiNegMom->SetParameters(325.6, 0.035, 0.175, 0.127, 0.139);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
            fParKaPosMom->SetParameters(15.3, 0.236, 0.203, 0.056, 0.494);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
            fParKaNegMom->SetParameters(8.88, 0.236, 0.203, 0.056, 0.494);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
            amplParam = 104.0;
            fParPrPosMom->SetParameters(amplParam, 0.213, 0.294, 0.09, 0.938);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
            amplParam /= fPrRatio;
            fParPrNegMom->SetParameters(amplParam, 0.213, 0.294, 0.09, 0.938);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
            fParDeMom->SetParameters(5.7, 0.338, 0.333, 0.114, 1.878);
            fPartYield[MpdPidUtils::kDeuteron].push_back(fParDeMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kDeuteron, fPartYield[MpdPidUtils::kDeuteron]));
            fParTrMom->SetParameters(0.2, -0.35, 0.723, 0.2, 2.81);
            fPartYield[MpdPidUtils::kTriton].push_back(fParTrMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kTriton, fPartYield[MpdPidUtils::kTriton]));
            fParHe3Mom->SetParameters(0.36, -0.784, 530.3, 0.131, 1.983);
            fPartYield[MpdPidUtils::kHe3].push_back(fParHe3Mom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe3, fPartYield[MpdPidUtils::kHe3]));
            fParHe4Mom->SetParameters(6.6e-03, 0.27, 0.2, 1.42, 3.51);
            fPartYield[MpdPidUtils::kHe4].push_back(fParHe4Mom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe4, fPartYield[MpdPidUtils::kHe4]));
         } else { /// fEnergy > 7.0 not ready, QGSM 9 gev
            fParElPosMom->SetParameters(17.6, -0.12, 0.078, 0.167, 0.00);
            fPartYield[MpdPidUtils::kElectron].push_back(fParElPosMom);
            fParElNegMom->SetParameters(16.3, -0.12, 0.078, 0.167, 0.00);
            fPartYield[MpdPidUtils::kElectron].push_back(fParElNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kElectron, fPartYield[MpdPidUtils::kElectron]));
            fParMuPosMom->SetParameters(20.5, 0.064, 0.107, 0.05, 0.105);
            fPartYield[MpdPidUtils::kMuon].push_back(fParMuPosMom);
            fParMuNegMom->SetParameters(20.5, 0.064, 0.107, 0.05, 0.105);
            fPartYield[MpdPidUtils::kMuon].push_back(fParMuNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kMuon, fPartYield[MpdPidUtils::kMuon]));
            fParPiPosMom->SetParameters(473., 0.034, 0.187, 0.469, 0.139);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
            fParPiNegMom->SetParameters(501.6, 0.034, 0.187, 0.469, 0.139);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
            fParKaPosMom->SetParameters(21.1, 0.157, 0.241, 0.043, 0.494);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
            fParKaNegMom->SetParameters(12.25, 0.157, 0.241, 0.043, 0.494);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
            amplParam = 67.4;
            fParPrPosMom->SetParameters(amplParam, .02, 0.365, 0.01, 0.938);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
            amplParam /= fPrRatio;
            fParPrNegMom->SetParameters(amplParam, .02, 0.365, 0.01, 0.938);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
            fParDeMom->SetParameters(1.8, 0.05, 0.432, 0.163, 1.878);
            fPartYield[MpdPidUtils::kDeuteron].push_back(fParDeMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kDeuteron, fPartYield[MpdPidUtils::kDeuteron]));
            fParTrMom->SetParameters(0.2, -0.35, 0.723, 0.2, 2.81);
            fPartYield[MpdPidUtils::kTriton].push_back(fParTrMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kTriton, fPartYield[MpdPidUtils::kTriton]));
            fParHe3Mom->SetParameters(0.36, -0.784, 530.3, 0.131, 1.983);
            fPartYield[MpdPidUtils::kHe3].push_back(fParHe3Mom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe3, fPartYield[MpdPidUtils::kHe3]));
            fParHe4Mom->SetParameters(6.6e-03, 0.27, 0.2, 1.42, 3.51);
            fPartYield[MpdPidUtils::kHe4].push_back(fParHe4Mom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe4, fPartYield[MpdPidUtils::kHe4]));
         }
      } else {                /// Tracking == "CF", "CFHM" is not ready
         if (fEnergy < 7.0) { /// not ready
            fParPiPosMom->SetParameters(44660.6, 1.4425, 0.0821871, 0.126819, -0.106069);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
            fParPiNegMom->SetParameters(1929.85, 24.5821, 0.111797, 0.066913, 0.667266);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
            fParKaPosMom->SetParameters(3089.81, 9.2121, 0.136135, 0.125175, 0.350274);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
            fParKaNegMom->SetParameters(237.81, -0.746975, 0.0771305, 0.0365433, 2.27772);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
            amplParam = 114392.;
            fParPrPosMom->SetParameters(amplParam, 6.35265, 0.188466, 0.166392, 0.7605);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
            amplParam /= fPrRatio;
            fParPrNegMom->SetParameters(amplParam, 6.35265, 0.188466, 0.166392, 0.7605);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
         } else { /// fEnergy = 11.0, 03/07/2018, 629K events
            fParPiPosMom->SetParameters(12708, 7.25584, 0.167878, 0.303888, 0.277294);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
            fParPiNegMom->SetParameters(3319.17, 14.7475, 0.212685, 0.341987, -0.190745);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
            fParKaPosMom->SetParameters(24333.9, 0.409004, 0.283049, 0.107653, 0.396297);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
            fParKaNegMom->SetParameters(12323.4, 0.392905, 0.267303, 0.0952731, 0.457674);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
            fParPrPosMom->SetParameters(42750.4, 1.92747, 0.363864, 0.273259, 0.337891);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
            fParPrNegMom->SetParameters(4636.8, 1.42461, 0.285596, 0.181618, 0.837135);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
         }
      }
   }
   if (Generator == "PHSD_CENT") { /// PHSD central events, b < 3 fm, fEnergy = 11.0, 20/03/2019, 5482 events
      fParPiPosMom->SetParameters(2628.47, 0.846834, 0.203337, 0.357357, 0.147232);
      fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
      fParPiNegMom->SetParameters(1997.0, 1.41096, 0.196863, 0.341347, 0.18051);
      fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
      fPartYieldMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
      fParKaPosMom->SetParameters(872.775, 0.346697, 0.298306, 0.102013, 0.423512);
      fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
      fParKaNegMom->SetParameters(437.025, 0.368078, 0.279169, 0.0955438, 0.48267);
      fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
      fPartYieldMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
      fParPrPosMom->SetParameters(1705.66, 0.663284, 0.425112, 0.220086, 0.457981);
      fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
      fParPrNegMom->SetParameters(191.997, 6.39986, 0.276029, 0.244932, 0.722006);
      fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
      fPartYieldMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
   }
   if (Generator == "PHSD_CSR") {
      if (fEnergy < 5.0) { /// PHSD csr central events, b < 3 fm, fEnergy = 4.0 GeV, 10/10/2019, ~50k events
         fParPiPosMom->SetParameters(326.591, 37.4694, 0.0395906, 0.145417, 1.08328);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
         fParPiNegMom->SetParameters(7834.26, 1.46868, 0.178981, 0.0634056, 0.520443);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
         fParKaPosMom->SetParameters(1709.8, 1.73908, 0.20194, 0.143219, 0.290164);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
         fParKaNegMom->SetParameters(327.052, 0.638412, 0.183259, 0.0800663, 0.554539);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
         fParPrPosMom->SetParameters(29176.4, 3.57695, 0.238379, 0.195002, 0.670283);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
         fParPrNegMom->SetParameters(4.18494, 1.07077, 0.165123, 0.0865973, 2.30082);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
      } else if (fEnergy < 7.0) { /// PHSD csr central events, b < 3 fm, fEnergy = 6.2 GeV, 10/10/2019, ~50k events
         fParPiPosMom->SetParameters(112.48, 188.057, 0.183905, 0.284109, -9.42658e-05);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
         fParPiNegMom->SetParameters(369.468, 63.2746, 0.190742, 0.278494, -0.0606429);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
         fParKaPosMom->SetParameters(4727.89, 0.411218, 0.271326, 0.10275, 0.38878);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
         fParKaNegMom->SetParameters(1467.24, 0.423877, 0.247845, 0.0908385, 0.477971);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
         fParPrPosMom->SetParameters(22653.3, 2.16885, 0.335439, 0.255945, 0.403967);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
         fParPrNegMom->SetParameters(85.4405, 1.87681, 0.13522, 0.0885255, 2.5824);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
      } else if (fEnergy < 8.0) { /// PHSD csr central events, b < 3 fm, fEnergy = 7.6 GeV, 10/10/2019, ~50k events
         fParPiPosMom->SetParameters(271.408, 92.8688, 0.184859, 0.303397, 0.0543418);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
         fParPiNegMom->SetParameters(1489.49, 17.1459, 0.192214, 0.302818, -0.0970675);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
         fParKaPosMom->SetParameters(6031.92, 0.356296, 0.28695, 0.0997535, 0.406052);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
         fParKaNegMom->SetParameters(2292.64, 0.392553, 0.26581, 0.0946355, 0.461694);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
         fParPrPosMom->SetParameters(20282.1, 1.38552, 0.373508, 0.255721, 0.370525);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
         fParPrNegMom->SetParameters(178.133, 1.51585, 0.201521, 0.125455, 1.68405);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
      } else if (fEnergy < 10.0) { /// PHSD csr central events, b < 3 fm, fEnergy = 8.8 GeV, 10/10/2019, ~50k events
         fParPiPosMom->SetParameters(330.227, 85.1517, 0.179918, 0.3141, 0.121203);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
         fParPiNegMom->SetParameters(21412.6, 0.965329, 0.210709, 0.338147, 0.08417);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
         fParKaPosMom->SetParameters(6714.74, 0.343555, 0.297204, 0.101954, 0.396783);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
         fParKaNegMom->SetParameters(2939.86, 0.373118, 0.276092, 0.0958796, 0.458351);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
         fParPrPosMom->SetParameters(18934.3, 1.02918, 0.397596, 0.248315, 0.378061);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
         fParPrNegMom->SetParameters(271.395, 1.29596, 0.246056, 0.146202, 1.31564);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
      } else { /// PHSD csr central events, b < 3 fm, fEnergy = 12.3 GeV, 10/10/2019, ~50k events
         fParPiPosMom->SetParameters(523.376, 66.8745, 0.172619, 0.330811, 0.198081);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
         fParPiNegMom->SetParameters(996.699, 25.9407, 0.202627, 0.35925, -0.162493);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
         fParKaPosMom->SetParameters(8208.95, 0.325369, 0.313731, 0.105367, 0.393535);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
         fParKaNegMom->SetParameters(4565.48, 0.344706, 0.298358, 0.100691, 0.441146);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
         fParPrPosMom->SetParameters(15861.8, 1.06114, 0.420781, 0.268906, 0.331336);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
         fParPrNegMom->SetParameters(545.455, 5.15299, 0.289769, 0.250636, 0.750215);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
      }
   }
   if (Generator == "PHSD_NOCSR") {
      if (fEnergy < 5.0) { /// PHSD no csr central events, b < 3 fm, fEnergy = 4 GeV, 10/10/2019, ~50k events
         fParPiPosMom->SetParameters(396.282, 32.7313, 0.0382901, 0.141021, 1.12647);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
         fParPiNegMom->SetParameters(3579.74, 3.98371, 0.11532, 0.0306925, 1.30859);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
         fParKaPosMom->SetParameters(1081.74, 3.09011, 0.190816, 0.154575, 0.279164);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
         fParKaNegMom->SetParameters(253.536, 0.589513, 0.178808, 0.0743182, 0.612162);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
         fParPrPosMom->SetParameters(30200.6, 3.29956, 0.236862, 0.190743, 0.684265);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
         fParPrNegMom->SetParameters(4.86563, 1.59102, 0.115168, 0.0729325, 3.3942);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
      } else if (fEnergy < 7.0) { /// PHSD no csr central events, b < 3 fm, fEnergy = 6.2 GeV, 10/10/2019, ~50k events
         fParPiPosMom->SetParameters(2907.1, 7.34039, 0.192883, 0.283338, 0.000931754);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
         fParPiNegMom->SetParameters(434.643, 57.0644, 0.194478, 0.275712, -0.0728979);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
         fParKaPosMom->SetParameters(3030.57, 0.431268, 0.264635, 0.102554, 0.388718);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
         fParKaNegMom->SetParameters(1194.21, 0.430582, 0.239526, 0.0874055, 0.502341);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
         fParPrPosMom->SetParameters(24464.4, 2.01624, 0.334706, 0.250743, 0.412746);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
         fParPrNegMom->SetParameters(104.891, 2.0149, 0.141495, 0.0947837, 2.53473);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
      } else if (fEnergy < 8.0) { /// PHSD no csr central events, b < 3 fm, fEnergy = 7.6 GeV, 10/10/2019, ~50k events
         fParPiPosMom->SetParameters(3463.66, 7.24569, 0.190893, 0.304263, 0.0839204);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
         fParPiNegMom->SetParameters(638.648, 40.2362, 0.194759, 0.301837, -0.1207);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
         fParKaPosMom->SetParameters(4243.36, 0.373601, 0.279679, 0.0997516, 0.401217);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
         fParKaNegMom->SetParameters(1954.45, 0.411682, 0.258866, 0.0941066, 0.470275);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
         fParPrPosMom->SetParameters(21992.6, 1.39463, 0.370946, 0.254295, 0.371529);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
         fParPrNegMom->SetParameters(210.959, 1.4898, 0.201572, 0.123791, 1.70536);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
      } else if (fEnergy < 10.0) { /// PHSD no csr central events, b < 3 fm, fEnergy = 8.8 GeV, 10/10/2019, ~50k events
         fParPiPosMom->SetParameters(6728.05, 3.88388, 0.190088, 0.319898, 0.124891);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
         fParPiNegMom->SetParameters(3648.45, 8.36849, 0.304356, 0.177689, 0.177577);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
         fParKaPosMom->SetParameters(5131.42, 0.356848, 0.289977, 0.101518, 0.39762);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
         fParKaNegMom->SetParameters(2595.77, 0.383127, 0.270511, 0.0949794, 0.464502);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
         fParPrPosMom->SetParameters(20165.6, 1.28612, 0.388469, 0.261465, 0.353915);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
         fParPrNegMom->SetParameters(305.183, 1.28427, 0.245721, 0.144331, 1.35441);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
      } else { /// PHSD no csr central events, b < 3 fm, fEnergy = 12.3 GeV, 10/10/2019, ~50k events
         fParPiPosMom->SetParameters(25403.9, 0.930567, 0.206834, 0.368397, 0.134282);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
         fParPiNegMom->SetParameters(1064.68, 36.2707, 0.326473, 0.173008, 0.228547);
         fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
         fParKaPosMom->SetParameters(7276.52, 0.328204, 0.309239, 0.104224, 0.392577);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
         fParKaNegMom->SetParameters(4294.11, 0.353445, 0.293973, 0.100527, 0.443493);
         fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
         fParPrPosMom->SetParameters(16438.3, 1.07511, 0.419069, 0.268969, 0.330482);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
         fParPrNegMom->SetParameters(575.422, 4.48476, 0.290547, 0.246113, 0.770269);
         fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
         fPartYieldMap.insert(
            pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
      }
   }
   if ((Generator == "LAQGSM") || (Generator == "QGSM")) {
      if (fTrackingState == MpdPidUtils::kHP) {
         if (fEnergy < 7.0) { /// not ready, QGSM 5 gev
            fParElPosMom->SetParameters(17.6, -0.12, 0.078, 0.167, 0.00);
            fPartYield[MpdPidUtils::kElectron].push_back(fParElPosMom);
            fParElNegMom->SetParameters(16.3, -0.12, 0.078, 0.167, 0.00);
            fPartYield[MpdPidUtils::kElectron].push_back(fParElNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kElectron, fPartYield[MpdPidUtils::kElectron]));
            fParMuPosMom->SetParameters(20.5, 0.064, 0.107, 0.05, 0.105);
            fPartYield[MpdPidUtils::kMuon].push_back(fParMuPosMom);
            fParMuNegMom->SetParameters(20.5, 0.064, 0.107, 0.05, 0.105);
            fPartYield[MpdPidUtils::kMuon].push_back(fParMuNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kMuon, fPartYield[MpdPidUtils::kMuon]));
            fParPiPosMom->SetParameters(307.0, 0.035, 0.175, 0.127, 0.139);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
            fParPiNegMom->SetParameters(325.6, 0.035, 0.175, 0.127, 0.139);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
            fParKaPosMom->SetParameters(15.3, 0.236, 0.203, 0.056, 0.494);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
            fParKaNegMom->SetParameters(8.88, 0.236, 0.203, 0.056, 0.494);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
            amplParam = 104.0;
            fParPrPosMom->SetParameters(amplParam, 0.213, 0.294, 0.09, 0.938);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
            amplParam /= fPrRatio;
            fParPrNegMom->SetParameters(amplParam, 0.213, 0.294, 0.09, 0.938);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
            fParDeMom->SetParameters(5.7, 0.338, 0.333, 0.114, 1.878);
            fPartYield[MpdPidUtils::kDeuteron].push_back(fParDeMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kDeuteron, fPartYield[MpdPidUtils::kDeuteron]));
            fParTrMom->SetParameters(0.2, -0.35, 0.723, 0.2, 2.81);
            fPartYield[MpdPidUtils::kTriton].push_back(fParTrMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kTriton, fPartYield[MpdPidUtils::kTriton]));
            fParHe3Mom->SetParameters(0.36, -0.784, 530.3, 0.131, 1.983);
            fPartYield[MpdPidUtils::kHe3].push_back(fParHe3Mom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe3, fPartYield[MpdPidUtils::kHe3]));
            fParHe4Mom->SetParameters(6.6e-03, 0.27, 0.2, 1.42, 3.51);
            fPartYield[MpdPidUtils::kHe4].push_back(fParHe4Mom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe4, fPartYield[MpdPidUtils::kHe4]));
         } else { /// fEnergy > 7.0 not ready, QGSM 9 gev
            fParElPosMom->SetParameters(17.6, -0.12, 0.078, 0.167, 0.00);
            fPartYield[MpdPidUtils::kElectron].push_back(fParElPosMom);
            fParElNegMom->SetParameters(16.3, -0.12, 0.078, 0.167, 0.00);
            fPartYield[MpdPidUtils::kElectron].push_back(fParElNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kElectron, fPartYield[MpdPidUtils::kElectron]));
            fParMuPosMom->SetParameters(20.5, 0.064, 0.107, 0.05, 0.105);
            fPartYield[MpdPidUtils::kMuon].push_back(fParMuPosMom);
            fParMuNegMom->SetParameters(20.5, 0.064, 0.107, 0.05, 0.105);
            fPartYield[MpdPidUtils::kMuon].push_back(fParMuNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kMuon, fPartYield[MpdPidUtils::kMuon]));
            fParPiPosMom->SetParameters(473., 0.034, 0.187, 0.469, 0.139);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
            fParPiNegMom->SetParameters(501.6, 0.034, 0.187, 0.469, 0.139);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
            fParKaPosMom->SetParameters(21.1, 0.157, 0.241, 0.043, 0.494);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
            fParKaNegMom->SetParameters(12.25, 0.157, 0.241, 0.043, 0.494);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
            amplParam = 67.4;
            fParPrPosMom->SetParameters(amplParam, .02, 0.365, 0.01, 0.938);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
            amplParam /= fPrRatio;
            fParPrNegMom->SetParameters(amplParam, .02, 0.365, 0.01, 0.938);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
            fParDeMom->SetParameters(1.8, 0.05, 0.432, 0.163, 1.878);
            fPartYield[MpdPidUtils::kDeuteron].push_back(fParDeMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kDeuteron, fPartYield[MpdPidUtils::kDeuteron]));
            fParTrMom->SetParameters(0.2, -0.35, 0.723, 0.2, 2.81);
            fPartYield[MpdPidUtils::kTriton].push_back(fParTrMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kTriton, fPartYield[MpdPidUtils::kTriton]));
            fParHe3Mom->SetParameters(0.36, -0.784, 530.3, 0.131, 1.983);
            fPartYield[MpdPidUtils::kHe3].push_back(fParHe3Mom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe3, fPartYield[MpdPidUtils::kHe3]));
            fParHe4Mom->SetParameters(6.6e-03, 0.27, 0.2, 1.42, 3.51);
            fPartYield[MpdPidUtils::kHe4].push_back(fParHe4Mom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe4, fPartYield[MpdPidUtils::kHe4]));
         }
      } else { /// Tracking == "CF", "CFHM" is not ready
         if (fEnergy < 7.0) {
            fParPiPosMom->SetParameters(44660.6, 1.4425, 0.0821871, 0.126819, -0.106069);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
            fParPiNegMom->SetParameters(1929.85, 24.5821, 0.111797, 0.066913, 0.667266);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
            fParKaPosMom->SetParameters(3089.81, 9.2121, 0.136135, 0.125175, 0.350274);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
            fParKaNegMom->SetParameters(237.81, -0.746975, 0.0771305, 0.0365433, 2.27772);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
            amplParam = 114392.;
            fParPrPosMom->SetParameters(amplParam, 6.35265, 0.188466, 0.166392, 0.7605);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
            amplParam /= fPrRatio;
            fParPrNegMom->SetParameters(amplParam, 6.35265, 0.188466, 0.166392, 0.7605);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
            fParDeMom->SetParameters(5024.67, 0.129733, 0.266767, 0.00308559, 39.0077);
            fPartYield[MpdPidUtils::kDeuteron].push_back(fParDeMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kDeuteron, fPartYield[MpdPidUtils::kDeuteron]));
            fParTrMom->SetParameters(938.334, 0.368862, 0.0161982, 0.00680544, 109.992);
            fPartYield[MpdPidUtils::kTriton].push_back(fParTrMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kTriton, fPartYield[MpdPidUtils::kTriton]));
            fParHe3Mom->SetParameters(661.212, 8.47325, 0.150273, 0.135588, 1.57518);
            fPartYield[MpdPidUtils::kHe3].push_back(fParHe3Mom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe3, fPartYield[MpdPidUtils::kHe3]));
         } else { /// fEnergy = 11.0 (mb, Au+Au), 12/12/2017, 394.5K events
            fParPiPosMom->SetParameters(1287.54, 48.9678, 0.253939, 0.121823, 0.371256);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
            fParPiNegMom->SetParameters(3979.89, 3.39065, 0.194877, 0.392421, -0.186646);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
            fParKaPosMom->SetParameters(9300.04, 0.252028, 0.280878, 0.0806716, 0.30038);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
            fParKaNegMom->SetParameters(4270.24, 0.269129, 0.276149, 0.0813025, 0.316033);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
            amplParam = 19730.6;
            fParPrPosMom->SetParameters(amplParam, 0.327478, 0.478771, 0.188208, 0.263743);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
            amplParam /= fPrRatio;
            fParPrNegMom->SetParameters(amplParam, 0.327478, 0.478771, 0.188208, 0.263743);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
            fParDeMom->SetParameters(303.123, 0.144, 0.580517, 0.097639, 1.31764);
            fPartYield[MpdPidUtils::kDeuteron].push_back(fParDeMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kDeuteron, fPartYield[MpdPidUtils::kDeuteron]));
         }
      }
   }
   if (Generator == "URQMD") {
      if (fTrackingState == MpdPidUtils::kHP) {
         if (fEnergy < 7.0) { /// not ready, QGSM 5 gev
            fParPiPosMom->SetParameters(307.0, 0.035, 0.175, 0.127, 0.139);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
            fParPiNegMom->SetParameters(325.6, 0.035, 0.175, 0.127, 0.139);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
            fParKaPosMom->SetParameters(15.3, 0.236, 0.203, 0.056, 0.494);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
            fParKaNegMom->SetParameters(8.88, 0.236, 0.203, 0.056, 0.494);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
            amplParam = 104.0;
            fParPrPosMom->SetParameters(amplParam, 0.213, 0.294, 0.09, 0.938);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
            amplParam /= fPrRatio;
            fParPrNegMom->SetParameters(amplParam, 0.213, 0.294, 0.09, 0.938);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
         } else { /// fEnergy > 7.0, fParameterization 11 GeV (27/07/2017) 10K events
            fPrRatio = 52.4698;
            fParPiPosMom->SetParameters(1050.82, 1.07492, 0.16913, 0.354773, 0.253396);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
            fParPiNegMom->SetParameters(1134.48, 1.05525, 0.170399, 0.35685, 0.241421);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
            fParKaPosMom->SetParameters(316.545, 0.193559, 0.359936, 0.0945873, 0.283851);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
            fParKaNegMom->SetParameters(174.667, 0.280725, 0.316105, 0.101935, 0.293752);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
            amplParam = 872.525;
            fParPrPosMom->SetParameters(amplParam, 0.440378, 0.475486, 0.207794, 0.440092);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
            amplParam /= fPrRatio;
            fParPrNegMom->SetParameters(amplParam, 0.440378, 0.475486, 0.207794, 0.440092);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
         }
      } else {
         if (fEnergy < 7.0) { /// not ready, QGSM 5 gev
            fParPiPosMom->SetParameters(307.0, 0.035, 0.175, 0.127, 0.139);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
            fParPiNegMom->SetParameters(325.6, 0.035, 0.175, 0.127, 0.139);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
            fParKaPosMom->SetParameters(15.3, 0.236, 0.203, 0.056, 0.494);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
            fParKaNegMom->SetParameters(8.88, 0.236, 0.203, 0.056, 0.494);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
            amplParam = 104.0;
            fParPrPosMom->SetParameters(amplParam, 0.213, 0.294, 0.09, 0.938);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
            amplParam /= fPrRatio;
            fParPrNegMom->SetParameters(amplParam, 0.213, 0.294, 0.09, 0.938);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
         } else { /// fEnergy > 7.0, fParameterization 8 GeV (27/07/2017) 10K events
            fPrRatio = 282.25;
            fParPiPosMom->SetParameters(3360.54, 1.08086, 0.151008, 0.330264, 0.351397);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
            fParPiNegMom->SetParameters(3530.8, 1.14233, 0.147187, 0.328013, 0.358974);
            fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
            fParKaPosMom->SetParameters(1206.47, 0.245233, 0.333437, 0.0996193, 0.327909);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
            fParKaNegMom->SetParameters(529.028, 0.288553, 0.303136, 0.0948728, 0.376275);
            fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
            amplParam = 4408.45;
            fParPrPosMom->SetParameters(amplParam, 0.666617, 0.416237, 0.215083, 0.487749);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
            amplParam /= fPrRatio;
            fParPrNegMom->SetParameters(amplParam, 0.666617, 0.416237, 0.215083, 0.487749);
            fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
            fPartYieldMap.insert(
               pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
         }
      }
   }
   if (Generator == "DEFAULT") { /// Generator == "DEFAULT", average 9 gev
      fParPiPosMom->SetParameters(685611., 3.51782, 0.285718, 0.187528, -0.141523);
      fPartYield[MpdPidUtils::kPion].push_back(fParPiPosMom);
      fParPiNegMom->SetParameters(427660., 8.83294, 0.281264, 0.194373, 0.0985359);
      fPartYield[MpdPidUtils::kPion].push_back(fParPiNegMom);
      fPartYieldMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kPion, fPartYield[MpdPidUtils::kPion]));
      fParKaPosMom->SetParameters(901126., 0.415003, 0.275805, 0.106005, 0.368604);
      fPartYield[MpdPidUtils::kKaon].push_back(fParKaPosMom);
      fParKaNegMom->SetParameters(398572., 0.409362, 0.261435, 0.0963171, 0.431606);
      fPartYield[MpdPidUtils::kKaon].push_back(fParKaNegMom);
      fPartYieldMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kKaon, fPartYield[MpdPidUtils::kKaon]));
      fParPrPosMom->SetParameters(1.92016e+06, 1.67276, 0.340098, 0.242819, 0.413844);
      fPartYield[MpdPidUtils::kProton].push_back(fParPrPosMom);
      fParPrNegMom->SetParameters(36391.7, 0.66405, 0.226386, 0.0955507, 1.41773);
      fPartYield[MpdPidUtils::kProton].push_back(fParPrNegMom);
      fPartYieldMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kProton, fPartYield[MpdPidUtils::kProton]));
      fParDeMom->SetParameters(38040.5, 0.129626, 0.078869, 0.0090554, 26.8717);
      fPartYield[MpdPidUtils::kDeuteron].push_back(fParDeMom);
      fPartYieldMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kDeuteron, fPartYield[MpdPidUtils::kDeuteron]));
      fParTrMom->SetParameters(855.693, 22.3505, 0.00224959, 0.00215819, 848.013);
      fPartYield[MpdPidUtils::kTriton].push_back(fParTrMom);
      fPartYieldMap.insert(
         pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kTriton, fPartYield[MpdPidUtils::kTriton]));
      fParHe3Mom->SetParameters(2705.9, -0.0458695, 0.469374, 0.455931, -0.456531);
      fPartYield[MpdPidUtils::kHe3].push_back(fParHe3Mom);
      fPartYieldMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe3, fPartYield[MpdPidUtils::kHe3]));
      fParHe4Mom->SetParameters(58.7835, 2.6726, 0.467089, 0.360178, 2.63697);
      fPartYield[MpdPidUtils::kHe4].push_back(fParHe4Mom);
      fPartYieldMap.insert(pair<MpdPidUtils::ePartType, vecTF1ptrs>(MpdPidUtils::kHe4, fPartYield[MpdPidUtils::kHe4]));
   }

   fGaus      = new TF1("fGaus", "gaus(0)", -1., 5.);
   fGaus2     = new TF2("fGaus2", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", -2., 10., -1., 5.);
   fAsymGaus  = new TF1("fAsymGaus", this, &MpdPid::AsymGaus, -1., 5., 4, "MpdPid", "AsymGaus");
   fAsymGaus2 = new TF2("fAsymGaus2", this, &MpdPid::AsymGaus2, -1., 5., -2., 10., 6, "MpdPid", "AsymGaus2");
}

vecTF1ptrs MpdPid::GetVecdEdxMean(MpdPidUtils::ePartType iType)
{
   vecTF1ptrs                                        ret;
   map<MpdPidUtils::ePartType, vecTF1ptrs>::iterator it = fdEdxBBMap.find(iType);
   if (it != fdEdxBBMap.end())
      ret = it->second;
   else
      cout << "Parameterization of <dE/dx> (" << MpdPidUtils::cParticleName[iType] << ") does not found." << endl;
   return ret;
}

vecTF1ptrs MpdPid::GetVecdEdxWidth(MpdPidUtils::ePartType iType)
{
   vecTF1ptrs                                        ret;
   map<MpdPidUtils::ePartType, vecTF1ptrs>::iterator it = fdEdxSigmaMap.find(iType);
   if (it != fdEdxSigmaMap.end())
      ret = it->second;
   else
      cout << "Parameterization of dE/dx width (" << MpdPidUtils::cParticleName[iType] << ") does not found." << endl;
   return ret;
}

vecTF1ptrs MpdPid::GetVecdEdxAsym(MpdPidUtils::ePartType iType)
{
   vecTF1ptrs                                        ret;
   map<MpdPidUtils::ePartType, vecTF1ptrs>::iterator it = fdEdxDeltaMap.find(iType);
   if (it != fdEdxDeltaMap.end())
      ret = it->second;
   else
      cout << "Parameterization of dE/dx delta parameter (" << MpdPidUtils::cParticleName[iType] << ") does not found."
           << endl;
   return ret;
}

vecTF1ptrs MpdPid::GetVecm2Width(MpdPidUtils::ePartType iType)
{
   vecTF1ptrs                                        ret;
   map<MpdPidUtils::ePartType, vecTF1ptrs>::iterator it = fParM2Map.find(iType);
   if (it != fParM2Map.end())
      ret = it->second;
   else
      cout << "Parameterization of m^2 width (" << MpdPidUtils::cParticleName[iType] << ") does not found." << endl;
   return ret;
}

vecTF1ptrs MpdPid::GetVecYield(MpdPidUtils::ePartType iType)
{
   vecTF1ptrs                                        ret;
   map<MpdPidUtils::ePartType, vecTF1ptrs>::iterator it = fPartYieldMap.find(iType);
   if (it != fPartYieldMap.end())
      ret = it->second;
   else
      cout << "Parameterization of yield (" << MpdPidUtils::cParticleName[iType] << ") does not found." << endl;
   return ret;
}

Double_t MpdPid::GetDedxParam(Double_t p, MpdPidUtils::ePartType iType)
{
   Double_t                                          dedx;
   map<MpdPidUtils::ePartType, vecTF1ptrs>::iterator ret = fdEdxBBMap.find(iType);
   if (ret != fdEdxBBMap.end()) {
      for (TF1 *BBF : ret->second) {
         if ((p < BBF->GetXmax()) && (p >= BBF->GetXmin())) {
            dedx = BBF->Eval(p);
            break;
         }
      }
   }

   return dedx;
}

Double_t MpdPid::Getm2WidthParam(Double_t p, MpdPidUtils::ePartType iType)
{
   Double_t                                          m2W;
   map<MpdPidUtils::ePartType, vecTF1ptrs>::iterator ret = fParM2Map.find(iType);
   if (ret != fParM2Map.end()) {
      for (TF1 *m2WF : ret->second) {
         if ((p < m2WF->GetXmax()) && (p >= m2WF->GetXmin())) {
            m2W = m2WF->Eval(p);
            break;
         }
      }
   }

   return m2W;
}

void MpdPid::ComputeBayesCoefficients(Double_t p)
{
   Int_t    maxType = fCharge == MpdPidUtils::kPos ? MpdPidUtils::kHe4 : MpdPidUtils::kProton;
   Double_t fsum    = 0.0;
   map<MpdPidUtils::ePartType, vecTF1ptrs>::iterator ret;

   for (Int_t iType = MpdPidUtils::kElectron; iType < MpdPidUtils::kNSpecies; iType++) fBayesCoefficients[iType] = 0.0;

   for (Int_t iType = MpdPidUtils::kElectron; iType < maxType + 1; iType++) {
      ret = fPartYieldMap.find(static_cast<MpdPidUtils::ePartType>(iType));
      if (ret != fPartYieldMap.end()) {
         fsum += ret->second[fCharge]->Eval(p);
      }
   }

   for (Int_t iType = MpdPidUtils::kElectron; iType < maxType + 1; iType++) {
      ret = fPartYieldMap.find(static_cast<MpdPidUtils::ePartType>(iType));
      if (ret != fPartYieldMap.end()) {
         fBayesCoefficients[iType] = ret->second[fCharge]->Eval(p) / fsum;
      }
   }
}

Bool_t MpdPid::FillProbs(MpdTrack *track)
{
   for (Int_t iType = 0; iType < MpdPidUtils::kNSpecies; iType++) fProb[iType] = 0.0;
   if (track == 0) return (kFALSE);
   Double_t px = track->GetPx(), py = track->GetPy(), pz = track->GetPz();
   Double_t p = TMath::Sqrt(px * px + py * py + pz * pz);
   ///
   /// NO parameterisation beyond p=5.0 GeV/c
   /// ignore anti-d,-He3 and -He4
   ///
   Int_t charge = 1;
   if (track->GetPt() > 0) charge = -1;
   Int_t    flag = track->GetTofFlag();
   Double_t dedx = track->GetdEdXTPC();
   Bool_t   ret;

   if (flag == 2 || flag == 6) {
      Double_t m2 = track->GetTofMass2();
      ret         = FillProbs(p, dedx, m2, charge);
   } else
      ret = FillProbs(p, dedx, charge);

   return ret;
}

Bool_t MpdPid::FillProbs(Double_t p, Double_t dedx, Int_t charge)
{
   ///
   /// NO parameterisation beyond p=5.0 GeV/c
   /// ignore anti-d,-He3 and -He4
   ///
   Double_t emean[MpdPidUtils::kNSpecies], sige[MpdPidUtils::kNSpecies], probs[MpdPidUtils::kNSpecies], fsum = 0.0, cut;

   for (Int_t iType = 0; iType < MpdPidUtils::kNSpecies; iType++) fProb[iType] = 0.0;

   if (p > 5.0) return kFALSE;
   if (charge == 0) return kFALSE;

   fCharge = charge > 0 ? MpdPidUtils::kPos : MpdPidUtils::kNeg;

   // if ( fSigmaEloss > 0.1 ) cut = fSigmaEloss;
   // else return kFALSE;

   Int_t nSigFlag = 0;

   ///
   /// Params. for dE/dx sigma versus total momentum.
   /// These are koef*<dE/dx> , koef=0.07
   /// Deviations from Bethe-Bloch at low-p accounted by somewhat larger sigmas
   ///

   for (Int_t iType = 0; iType < MpdPidUtils::kNSpecies; iType++) {
      emean[iType] = GetDedxParam(p, static_cast<MpdPidUtils::ePartType>(iType));
      sige[iType]  = GetDedxWidthValue(p, static_cast<MpdPidUtils::ePartType>(iType));
   }

   if (fMethod) { /// bayesian approach
      ComputeBayesCoefficients(p);
      ///
      /// Set prob=0. for differences greater than "n" of sigmas
      /// otherwise evaluation..
      ///

      for (Int_t iType = 0; iType < MpdPidUtils::kNSpecies; iType++) {
         if (dedx < emean[iType])
            cut = fdEdxBottomBound[iType];
         else
            cut = fdEdxTopBound[iType];
         if (cut > 0.1)
            probs[iType] = ComputeDedxProb_asym(cut, p, dedx, fBayesCoefficients[iType], emean[iType], sige[iType],
                                                static_cast<MpdPidUtils::ePartType>(iType));
         else
            probs[iType] = 0.0;
         fsum += probs[iType];
      }

      ///
      /// Normalization
      ///

      if (1.0e+8 * fsum > 0.0) {
         for (Int_t iType = 0; iType < MpdPidUtils::kNSpecies; iType++) fProb[iType] = probs[iType] / fsum;
      } else
         return kFALSE; /// outliers!

      return kTRUE;
   } else { /// n-sigma method
      for (Int_t iType = 0; iType < MpdPidUtils::kNSpecies; iType++) {
         ComputeEnLossSigma(p, dedx, emean[iType], sige[iType], static_cast<MpdPidUtils::ePartType>(iType));
         if (dedx < emean[iType])
            cut = fdEdxBottomBound[iType];
         else
            cut = fdEdxTopBound[iType];
         if ((fEnLossSigmasArray[iType] < cut) && (fNSigSpecies[iType])) {
            fProb[iType] = 1.0;
            nSigFlag++;
         } else
            fProb[iType] = 0.0;
      }
      if (nSigFlag > 1) {
         for (Int_t iType = 0; iType < MpdPidUtils::kNSpecies; iType++) fProb[iType] = 0.0;
      }
      if (nSigFlag == 1)
         return kTRUE;
      else
         return kFALSE;
   }
}

Bool_t MpdPid::FillProbs(Double_t p, Double_t dedx, Double_t m2, Int_t charge)
{
   ///
   /// NO parameterisation beyond p=3.0 GeV/c
   /// ignore anti-d,-He3 and -He4
   ///
   Double_t emean[MpdPidUtils::kNSpecies], sige[MpdPidUtils::kNSpecies], sigm[MpdPidUtils::kNSpecies],
      probs[MpdPidUtils::kNSpecies], fsum = 0.0, cut_dedx, cut_m2;

   for (Int_t iType = 0; iType < MpdPidUtils::kNSpecies; iType++) fProb[iType] = 0.0;

   if (p > 5.0) return kFALSE;
   if (charge == 0) return kFALSE;

   // if ( (fSigmaEloss > 0.1) && (fSigmaTof > 0.1) )
   //{ cut_dedx = fSigmaEloss; cut_m2 = fSigmaTof; }
   // else return kFALSE;

   fCharge = charge > 0 ? MpdPidUtils::kPos : MpdPidUtils::kNeg;

   Int_t nSigFlag = 0;

   ///
   /// Params. for dE/dx sigma versus total momentum.
   /// These are koef*<dE/dx> , koef=0.06
   /// Deviations from Bethe-Bloch at low-p accounted by somewhat larger sigmas
   ///

   for (Int_t iType = 0; iType < MpdPidUtils::kNSpecies; iType++) {
      emean[iType] = GetDedxParam(p, static_cast<MpdPidUtils::ePartType>(iType));
      sige[iType]  = GetDedxWidthValue(p, static_cast<MpdPidUtils::ePartType>(iType));
      sigm[iType]  = Getm2WidthParam(p, static_cast<MpdPidUtils::ePartType>(iType));
   }

   ///
   /// Set prob=0. for differences greater than n-sigmas
   /// otherwise evaluation..
   ///

   const Double_t aMass[MpdPidUtils::kNSpecies] = {0.0007, 0.011, 0.019, 0.24, 0.887, 3.54, 7.87, 1.983, 3.51};

   if (fMethod) { /// bayesian approach
      ComputeBayesCoefficients(p);
      for (Int_t iType = 0; iType < MpdPidUtils::kNSpecies; iType++) {
         if (dedx < emean[iType])
            cut_dedx = fdEdxBottomBound[iType];
         else
            cut_dedx = fdEdxTopBound[iType];
         if (m2 < aMass[iType])
            cut_m2 = fm2BottomBound[iType];
         else
            cut_m2 = fm2TopBound[iType];
         probs[iType] =
            ComputeCombProb_asym(cut_dedx, cut_m2, p, dedx, m2, fBayesCoefficients[iType], emean[iType], aMass[iType],
                                 sige[iType], sigm[iType], static_cast<MpdPidUtils::ePartType>(iType));
         fsum += probs[iType];
      }

      ///
      /// Normalization
      ///

      if (1.0e+8 * fsum > 0.0) {
         for (Int_t iType = 0; iType < MpdPidUtils::kNSpecies; iType++) fProb[iType] = probs[iType] / fsum;
      } else
         return kFALSE; /// outliers!

      return kTRUE;
   } else { /// n-sigma method
      for (Int_t iType = 0; iType < MpdPidUtils::kNSpecies; iType++) {
         ComputeEnLossSigma(p, dedx, emean[iType], sige[iType], static_cast<MpdPidUtils::ePartType>(iType));
         ComputeMSquaredSigma(m2, aMass[iType], sigm[iType], static_cast<MpdPidUtils::ePartType>(iType));
         if (dedx < emean[iType])
            cut_dedx = fdEdxBottomBound[iType];
         else
            cut_dedx = fdEdxTopBound[iType];
         if (m2 < aMass[iType])
            cut_m2 = fm2BottomBound[iType];
         else
            cut_m2 = fm2TopBound[iType];

         // if ( ( TMath::Sqrt( TMath::Power( fEnLossSigmasArray[iType], 2) + TMath::Power( fMSquaredSigmasArray[iType],
         // 2) ) <= TMath::Sqrt( TMath::Power( cut_dedx, 2 ) + TMath::Power( cut_m2, 2 )) ) && (fNSigSpecies[iType]) ) {
         if ((fEnLossSigmasArray[iType] <= cut_dedx) && (fMSquaredSigmasArray[iType] < cut_m2) &&
             (fNSigSpecies[iType])) {
            fProb[iType] = 1.0;
            nSigFlag++;
         } else
            fProb[iType] = 0.0;
      }
      if (nSigFlag > 1) {
         for (Int_t iType = 0; iType < MpdPidUtils::kNSpecies; iType++) fProb[iType] = 0.0;
      }
      if (nSigFlag == 1)
         return kTRUE;
      else
         return kFALSE;
   }
}

Long_t MpdPid::GetMaxProb()
{
   Long_t   pdg                           = 211;
   Double_t pcut                          = 0.501;
   Long_t   codes[MpdPidUtils::kNSpecies] = {-11, -13, 211, 321, 2212, 1000010020, 1000010030, 1000020030, 1000020040};
   Double_t probs[MpdPidUtils::kNSpecies] = {
      fProb[MpdPidUtils::kElectron], fProb[MpdPidUtils::kMuon],   fProb[MpdPidUtils::kPion],
      fProb[MpdPidUtils::kKaon],     fProb[MpdPidUtils::kProton], fProb[MpdPidUtils::kDeuteron],
      fProb[MpdPidUtils::kTriton],   fProb[MpdPidUtils::kHe3],    fProb[MpdPidUtils::kHe4]};

   if (fProb[MpdPidUtils::kElectron] > pcut)
      pdg = -11;
   else if (fProb[MpdPidUtils::kMuon] > pcut)
      pdg = -13;
   else if (fProb[MpdPidUtils::kPion] > pcut)
      pdg = 211;
   else if (fProb[MpdPidUtils::kKaon] > pcut)
      pdg = 321;
   else if (fProb[MpdPidUtils::kProton] > pcut)
      pdg = 2212;
   else if (fProb[MpdPidUtils::kDeuteron] > pcut)
      pdg = 1000010020;
   else if (fProb[MpdPidUtils::kTriton] > pcut)
      pdg = 1000010030;
   else if (fProb[MpdPidUtils::kHe3] > pcut)
      pdg = 1000020030;
   else if (fProb[MpdPidUtils::kHe4] > pcut)
      pdg = 1000020040;
   else {
      Long_t tmp;
      for (Int_t i = 1; i < MpdPidUtils::kNSpecies; i++) {
         if (probs[0] < probs[i]) {
            tmp      = codes[0];
            codes[0] = codes[i];
            codes[i] = tmp;
            tmp      = probs[0];
            probs[0] = probs[i];
            probs[i] = tmp;
         }
      }
      pdg = codes[0];
   }
   pdg *= (fCharge == MpdPidUtils::kPos) ? 1.0 : -1.0;

   return pdg;
}

Double_t MpdPid::GetNsigmaToBetheBloch(TString species)
{
   Double_t ans = -1.0;
   if (species.EqualTo("el"))
      ans = fEnLossSigmasArray[MpdPidUtils::kElectron];
   else if (species.EqualTo("mu"))
      ans = fEnLossSigmasArray[MpdPidUtils::kMuon];
   else if (species.EqualTo("pi"))
      ans = fEnLossSigmasArray[MpdPidUtils::kPion];
   else if (species.EqualTo("ka"))
      ans = fEnLossSigmasArray[MpdPidUtils::kKaon];
   else if (species.EqualTo("pr"))
      ans = fEnLossSigmasArray[MpdPidUtils::kProton];
   else if (species.EqualTo("de"))
      ans = fEnLossSigmasArray[MpdPidUtils::kDeuteron];
   else if (species.EqualTo("tr"))
      ans = fEnLossSigmasArray[MpdPidUtils::kTriton];
   else if (species.EqualTo("he3"))
      ans = fEnLossSigmasArray[MpdPidUtils::kHe3];
   else if (species.EqualTo("he4"))
      ans = fEnLossSigmasArray[MpdPidUtils::kHe4];
   return ans;
}

Double_t MpdPid::GetNsigmaToBetheBloch(MpdPidUtils::ePartType iType)
{
   Double_t ans = (iType != MpdPidUtils::kUnknown) ? fEnLossSigmasArray[iType] : -1.0;
   return ans;
}

Double_t MpdPid::GetNsigmaToAverageMass2(TString species)
{
   Double_t ans = -1.0;
   if (species.EqualTo("el"))
      ans = fMSquaredSigmasArray[MpdPidUtils::kElectron];
   else if (species.EqualTo("mu"))
      ans = fMSquaredSigmasArray[MpdPidUtils::kMuon];
   else if (species.EqualTo("pi"))
      ans = fMSquaredSigmasArray[MpdPidUtils::kPion];
   else if (species.EqualTo("ka"))
      ans = fMSquaredSigmasArray[MpdPidUtils::kKaon];
   else if (species.EqualTo("pr"))
      ans = fMSquaredSigmasArray[MpdPidUtils::kProton];
   else if (species.EqualTo("de"))
      ans = fMSquaredSigmasArray[MpdPidUtils::kDeuteron];
   else if (species.EqualTo("tr"))
      ans = fMSquaredSigmasArray[MpdPidUtils::kTriton];
   else if (species.EqualTo("he3"))
      ans = fMSquaredSigmasArray[MpdPidUtils::kHe3];
   else if (species.EqualTo("he4"))
      ans = fMSquaredSigmasArray[MpdPidUtils::kHe4];
   return ans;
}

Double_t MpdPid::GetNsigmaToAverageMass2(MpdPidUtils::ePartType iType)
{
   Double_t ans = (iType != MpdPidUtils::kUnknown) ? fMSquaredSigmasArray[iType] : -1.0;
   return ans;
}

void MpdPid::SetPrRat(Double_t PrRat)
{
   fPrRatio                                             = PrRat;
   map<MpdPidUtils::ePartType, vecTF1ptrs>::iterator it = fPartYieldMap.find(MpdPidUtils::kProton);
   if (it != fPartYieldMap.end()) {
      it->second[1]->SetParameter(0, it->second[1]->GetParameter(0) / PrRat);
   }
}

Double_t MpdPid::GetBayesCoefficient(Double_t p, MpdPidUtils::ePartType iType, MpdPidUtils::ePartCharge ech)
{
   Int_t    maxType = ech == MpdPidUtils::kPos ? MpdPidUtils::kHe4 : MpdPidUtils::kProton;
   Double_t fsum    = 0.0;
   map<MpdPidUtils::ePartType, vecTF1ptrs>::iterator ret;

   for (Int_t itype = MpdPidUtils::kElectron; itype < MpdPidUtils::kNSpecies; itype++) fBayesCoefficients[itype] = 0.0;

   for (Int_t itype = MpdPidUtils::kElectron; itype < maxType + 1; itype++) {
      ret = fPartYieldMap.find(static_cast<MpdPidUtils::ePartType>(itype));
      if (ret != fPartYieldMap.end()) {
         fsum += ret->second[ech]->Eval(p);
      }
   }

   Double_t ans = -1.0;
   ret          = fPartYieldMap.find(static_cast<MpdPidUtils::ePartType>(iType));
   if (ret != fPartYieldMap.end())
      ans = ret->second[ech]->Eval(p) / fsum;
   else
      cout << "Bayes Coefficient cannot be computed." << endl;
   return ans;
}

void MpdPid::Print(const char *comment, ostream &os)
{
   if (comment != nullptr) os << comment;

   os << "MpdPid::Print" << endl;
   for (Int_t iType = MpdPidUtils::kElectron; iType < MpdPidUtils::kNSpecies + 1; iType++)
      os << "fProb[" << MpdPidUtils::cParticleShortName[iType] << "] = " << fProb[iType] << endl;
}

ClassImp(MpdPid);
