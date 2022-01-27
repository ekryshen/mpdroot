
//  Class that holds TPC Hit information.
//
//  Original author: Jonathan Paley   08/18/04   jpaley@indiana.edu
//  Edited by Artem Basalaev 03/18/13

#ifndef MPDTPCFOUNDHIT_H
#define MPDTPCFOUNDHIT_H

#include <iostream>
#include <iomanip>
#include "TObject.h"
#include "TRef.h"
#include "MpdTpc2dCluster.h"

using namespace std;

class MpdTpcFoundHit : public TObject {
public:
   typedef enum { kUnknown = 0, kWMPeak, kFitPeak, kMixFitWM } TpcFoundHitType_t;

public:
   MpdTpcFoundHit();
   MpdTpcFoundHit(Float_t x, Float_t dx, Float_t y, Float_t dy, Float_t z, Float_t dz);
   MpdTpcFoundHit(Float_t p[3], Float_t dp[3]);
   MpdTpcFoundHit(const MpdTpcFoundHit &hit);
   MpdTpcFoundHit(const MpdTpc2dCluster *clus);
   ~MpdTpcFoundHit();

public:
   const MpdTpc2dCluster *Cluster() const
   {
      // return dynamic_cast<const MpdTpc2dCluster*>(fCluster.GetObject());
      return fCluster;
   }
   void AttachCluster(const MpdTpc2dCluster *);

   const Int_t   PadRow() const { return fPadRow; }
   const Float_t PadCol() const { return fPadCol; }
   const Float_t TimeBkt() const { return fTimeBkt; }

   const Float_t QFit() const { return fQFit; }
   const Float_t QFitSig() const { return fSigQFit; }
   const Float_t QADC() const { return fQADC; }

   const Float_t GetLocalX() const { return _xl; }
   const Float_t GetLocalY() const { return _yl; }
   const Float_t GetLocalZ() const { return _zl; }
   const Float_t GetGlobalX() const { return _xg; }
   const Float_t GetGlobalY() const { return _yg; }
   const Float_t GetGlobalZ() const { return _zg; }
   const Float_t errX() const { return _dx; }
   const Float_t errY() const { return _dy; }
   const Float_t errZ() const { return _dz; }

   const Int_t MaxCol() const { return this->Cluster()->MaxCol(); }
   const Int_t MinCol() const { return this->Cluster()->MinCol(); }
   const Int_t MaxBkt() const { return this->Cluster()->MaxBkt(); }
   const Int_t MinBkt() const { return this->Cluster()->MinBkt(); }
   const Int_t GetSect() const { return this->Cluster()->GetSect(); }

   const Int_t NumHits() const { return fNumHits; }
   const Int_t Type() const { return fHitType; }

   const Int_t GetOrigin() const { return fOrigin; }

   void SetTimeBkt(Float_t val) { fTimeBkt = val; }
   void SetPadCol(Float_t col) { fPadCol = col; }
   void SetNumHits(Int_t val) { fNumHits = val; }
   void SetType(Int_t val) { fHitType = val; }
   void SetQFit(Float_t val) { fQFit = val; }
   void SetSigQFit(Float_t val) { fSigQFit = val; }
   void SetQADC(Float_t val) { fQADC = val; }
   void SetPos(Float_t *p, Float_t *dp);
   void SetGlobalX(Float_t x) { _xg = x; }
   void SetGlobalY(Float_t y) { _yg = y; }
   void SetGlobalZ(Float_t z) { _zg = z; }
   void SetXerr(Float_t errx) { _dx = errx; }
   void SetYerr(Float_t erry) { _dy = erry; }
   void SetZerr(Float_t errz) { _dz = errz; }
   void AddOrigin(Int_t Origin);

   const void Print() const;

public:
   friend ostream &operator<<(ostream &ostr, const MpdTpcFoundHit &h1);

private:
   Int_t   fOrigin;
   short   fNumHits;
   short   fHitType;
   short   fPadRow;  // Pad row; converted to y-position
   Float_t fPadCol;  // Pad column; converted to x-position
   Float_t fTimeBkt; // Time bucket;
   // QFit is the total charge from cluster slices fit to Gamma
   // (area under fit curve)
   Float_t fQFit;
   // SigQFit is the "uncertainty" on the total charge from the Gamma fit
   // Note: this isn't a real uncertainty, since a true chi^2 was not calculated
   // but rather a goodness-of-fit
   Float_t fSigQFit;
   // QADC is the total charge from cluster slices that we failed to fit to
   // a Gamma (sum of ADC values)
   Float_t fQADC;

   // TRef fCluster; // poInt_ter to the cluster from which hit was formed
   const MpdTpc2dCluster *fCluster; // poInt_ter to the cluster from which hit was formed
   Float_t                _xg, _yg, _zg, _xl, _yl, _zl, _dx, _dy, _dz;

   ClassDef(MpdTpcFoundHit, 4)
};

#endif // TPCRHIT
