#ifndef LHE_HIT_H
#define LHE_HIT_H

// *************************************************************************
//  Author: Oleg Rogachevsky             e-mail: rogach@sunhe.jinr.ru
//
// hit information for LHE
//
// Created: 1-07-07
// Modified:
//
// *************************************************************************

#include "TObject.h"
#include "TClonesArray.h"
#include "TVector3.h"

#include "TpcPoint.h"

class TpcLheTrack;

class TpcLheHit : public TObject {

private:
   TVector3 fCoord;     // vector of hit coordinates
   TVector3 fError;     // vector of errors on hit coordinates
   Int_t    fHitNumber; // number of this hit in this event
   Int_t    fTrackID;   // track number to which this hit belongs to
   Int_t    fRefIndex;  // ref index on MC point
   Bool_t   fUsed;      // if the hit is assigned to a track
   Double_t fR;         // hit radius
   Double_t fRphi;      // hit R-Phi coordinate
   Int_t    fFlag;      // flag
   Int_t    fLayer;     // TPC layer
   Double_t fEdep;      // energy deposit
   Double_t fStep;      // step size

public:
   TpcLheHit();                       // default constructor
   TpcLheHit(TpcPoint *point);        // constructor for data after hit finding
   TpcLheHit(Double_t *x, Int_t stn); // constructor which take an arbitrary point
   virtual ~TpcLheHit();              // destructor

   // getters
   TVector3 GetCoord() { return fCoord; }
   TVector3 GetError() { return fError; }

   Double_t GetX() const { return fCoord.X(); }
   Double_t GetY() const { return fCoord.Y(); }
   Double_t GetZ() const { return fCoord.Z(); }
   Double_t GetXerr() const { return fError.X(); }
   Double_t GetYerr() const { return fError.Y(); }
   Double_t GetZerr() const { return fError.Z(); }
   Double_t GetRphiErr() const { return 0.05; } // AZ 0.5 mm error in RPhi
   Double_t GetR() const { return fR; }
   Double_t GetRphi() const { return fRphi; }
   Int_t    GetFlag() const { return fFlag; }
   Int_t    GetLayer() const { return fLayer; }
   Double_t GetEdep() const { return fEdep; }
   Double_t GetStep() const { return fStep; }

   // TpcLheTrack *GetTrack(TClonesArray *tracks) const;

   Bool_t GetUsage() const { return fUsed; }
   Int_t  GetHitNumber() const { return fHitNumber; }
   Int_t  GetTrackID() const { return fTrackID; }
   Int_t  GetRefIndex() const { return fRefIndex; }

   // setters
   void SetX(Double_t f) { fCoord.SetX(f); }
   void SetY(Double_t f) { fCoord.SetY(f); }
   void SetZ(Double_t f) { fCoord.SetZ(f); }
   void SetXerr(Double_t f) { fError.SetX(f); }
   void SetYerr(Double_t f) { fError.SetY(f); }
   void SetZerr(Double_t f) { fError.SetZ(f); }
   void SetR(Double_t f) { fR = f; }
   void SetRphi(Double_t f) { fRphi = f; }
   void SetFlag(Int_t f) { fFlag = f; }
   void SetLayer(Int_t f) { fLayer = f; }
   void SetEdep(Double_t edep) { fEdep = edep; }
   void SetStep(Double_t step) { fStep = step; }

   void SetUsage(Bool_t f) { fUsed = f; }
   void SetHitNumber(Int_t f) { fHitNumber = f; }
   void SetTrackID(Int_t f) { fTrackID = f; }
   void SetRefIndex(Int_t k) { fRefIndex = k; }

   void Print();

   Bool_t IsSortable() const { return kTRUE; }
   Int_t  Compare(const TObject *hit) const; // "Compare" function for sorting

   ClassDef(TpcLheHit, 1) //
};

#endif
