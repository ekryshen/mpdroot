#include "MpdTpcPeak.h"
#include <cmath>

ClassImp(MpdTpcPeak);

MpdTpcPeak::MpdTpcPeak()
   : fCol(-1), fBktOff(0), fPeakTime(0.), fIntegral(0.), fIntegSig(0.), fChi2(-1.), fMean(-1.e5), fSigMean(-1.),
     fSumADC(-1.), fMaxADC(0.), fMinADC(1.e9)
{
   fADC.clear();
}

//......................................................................

MpdTpcPeak::MpdTpcPeak(Int_t col, Float_t ptime, Float_t integ, Float_t chi2)
{
   fCol      = col;
   fPeakTime = ptime;
   fIntegral = integ;
   fChi2     = chi2;
   fBktOff   = 0;
   fMean     = -1.;
   fSigMean  = -1.;
   fSumADC   = -1.;
   fMaxADC   = 0.;
   fMinADC   = 1.e9;
   fADC.clear();
}

//......................................................................

MpdTpcPeak::MpdTpcPeak(const MpdTpcPeak &peak) : TObject((const TObject &)peak)
{

   fCol      = peak.fCol;
   fBktOff   = peak.fBktOff;
   fCluster  = peak.fCluster;
   fChi2     = peak.fChi2;
   fPeakTime = peak.fPeakTime;
   fIntegral = peak.fIntegral;
   fIntegSig = peak.fIntegSig;
   fMean     = peak.fMean;
   fSigMean  = peak.fSigMean;
   fMaxADC   = peak.fMaxADC;
   fMinADC   = peak.fMinADC;
   fSumADC   = 0.;

   for (UInt_t i = 0; i < peak.fADC.size(); ++i) {
      fADC.push_back(peak.fADC[i]);
      fSumADC += peak.fADC[i];
   }
}

//......................................................................

MpdTpcPeak::~MpdTpcPeak()
{
   fADC.clear();
}

//......................................................................

void MpdTpcPeak::Insert(Float_t val, Int_t col, Int_t bkt)
{
   fADC.push_back(val);
   if (val > fMaxADC) {
      fMaxADC    = val;
      fMaxAdcCol = col;
      fMaxAdcBkt = bkt;
   }
   if (val < fMinADC) fMinADC = val;
}

//......................................................................

void MpdTpcPeak::Remove(UInt_t val)
{
   fADC.erase(fADC.begin() + val);
   fMinADC = 1.e9;
   fMaxADC = 0.;

   for (UInt_t i = 0; i < fADC.size(); ++i) {
      if (fADC[i] > fMaxADC) fMaxADC = fADC[i];
      if (fADC[i] < fMinADC) fMinADC = fADC[i];
   }
}

//......................................................................

void MpdTpcPeak::Clear()
{
   fCol      = -1;
   fCluster  = 0;
   fChi2     = -1.;
   fPeakTime = 0.;
   fIntegral = 0.;
   fIntegSig = 0.;
   fBktOff   = 0;
   fMean     = -1.e5;
   fSigMean  = -1.;
   fSumADC   = -1.;
   fMaxADC   = 0.;
   fMinADC   = 1.e9;

   fADC.clear();
}

//......................................................................

void MpdTpcPeak::AttachCluster(const MpdTpc2dCluster *currclus)
{
   // fCluster = TRef((MpdTpc2dCluster*)currclus);
   fCluster = currclus;
}

//......................................................................

const Float_t MpdTpcPeak::Mean()
{

   if (fMean < -1.e2) {
      fMean        = 0.;
      Float_t sumq = 0.;
      for (UInt_t i = 0; i < fADC.size(); ++i) {
         sumq += fADC[i];
         fMean += fADC[i] * (Float_t)i;
      }
      fMean /= sumq;
   }

   return fMean;
}

//......................................................................

const Float_t MpdTpcPeak::SigMean()
{

   if (fSigMean < 0.) {
      Float_t mean = this->Mean();
      fSigMean     = 0.;
      for (UInt_t i = 0; i < fADC.size(); ++i) {
         fSigMean += (fADC[i] - mean) * (fADC[i] - mean);
      }
      fSigMean /= (Float_t)fADC.size();
      fSigMean = sqrt(fSigMean);
      fSigMean /= sqrt((Float_t)fADC.size() - 1);
   }

   return fSigMean;
}

//......................................................................

const Float_t MpdTpcPeak::SumADC()
{

   if (fSumADC < 0.) {
      fSumADC = 0.;
      for (UInt_t i = 0; i < fADC.size(); ++i) fSumADC += fADC[i];
   }

   return fSumADC;
}

////////////////////////////////////////////////////////////////////////
