
//  Class that holds TPC Pad Peak information.
//
//  Original author: Jonathan Paley   05/13/06   jpaley@indiana.edu
//  Edited by Artem Basalaev 03/18/13

#ifndef MPDTPCPEAK_H
#define MPDTPCPEAK_H
#include <iostream>
#include <iomanip>
#include "TObject.h"
#include "TRef.h"
#include "MpdTpc2dCluster.h"

class MpdTpcPeak : public TObject {
public:
   MpdTpcPeak();
   MpdTpcPeak(Int_t col, Float_t ptime, Float_t integ, Float_t chi2);
   MpdTpcPeak(const MpdTpcPeak &peak);
   ~MpdTpcPeak();

public:
   const MpdTpc2dCluster *Cluster() const
   {
      // return dynamic_cast<const MpdTpc2dCluster*>(fCluster.GetObject());
      return fCluster;
   }
   void AttachCluster(const MpdTpc2dCluster *);

   const Int_t   Row() const { return this->Cluster()->Row(); }
   const Int_t   Col() const { return fCol; }
   const Int_t   BktOff() const { return fBktOff; }
   const Int_t   NSamples() const { return (Int_t)fADC.size(); }
   const Float_t Sample(Int_t i) const { return fADC[i]; }
   const Float_t PeakTime() const { return fPeakTime; }
   const Float_t Integral() const { return fIntegral; }
   const Float_t IntegSig() const { return fIntegSig; }
   const Float_t Chi2() const { return fChi2; }
   const Float_t Max() const { return fMaxADC; }
   const Int_t   MaxAdcCol() const { return fMaxAdcCol; }
   const Int_t   MaxAdcBkt() const { return fMaxAdcBkt; }
   const Float_t Min() const { return fMinADC; }
   const Int_t   GetOrigin() const { return fOrigin; }
   const Float_t Mean();
   const Float_t SigMean();
   const Float_t SumADC();

   void Insert(Float_t val, Int_t col, Int_t bkt);
   void Remove(UInt_t val);
   void SetCol(Int_t val) { fCol = val; }
   void SetBktOff(Int_t val) { fBktOff = val; }
   void SetPeakTime(Float_t val) { fPeakTime = val; }
   void SetIntegral(Float_t val) { fIntegral = val; }
   void SetIntegSig(Float_t val) { fIntegSig = val; }
   void SetChi2(Float_t val) { fChi2 = val; }
   void SetOrigin(Int_t Origin) { fOrigin = Origin; }
   void Clear();

private:
   Int_t                fCol;
   Int_t                fBktOff;
   std::vector<Float_t> fADC;
   Float_t              fPeakTime;
   Float_t              fIntegral;
   Float_t              fIntegSig;
   Float_t              fChi2;
   Float_t              fMean;
   Float_t              fSigMean;
   Float_t              fSumADC;
   Float_t              fMaxADC;
   Float_t              fMinADC;
   Int_t                fMaxAdcCol;
   Int_t                fMaxAdcBkt;
   Int_t                fOrigin;

   // TRef fCluster; // pointer to the cluster from which hit was formed
   const MpdTpc2dCluster *fCluster; // pointer to the cluster from which hit was formed

   ClassDef(MpdTpcPeak, 1)
};

#endif // TPCRPADPEAK
