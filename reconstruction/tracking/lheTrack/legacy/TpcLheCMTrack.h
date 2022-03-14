/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// TpcLheCMTrack class - representation of one STS track with CM points        //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////

#ifndef PND_STS_CM_Track_H
#define PND_STS_CM_Track_H

#include "TObject.h"
#include "TObjArray.h"
#include "TVector3.h"

#include "TpcLheTrack.h"
#include "TpcLheCMPoint.h"
#include "TpcLheHit.h"

#include "Riostream.h"

class TpcLheCMTrack : public TpcLheTrack {

private:
   Int_t      fNMapHits;
   TObjArray *fMappingHits; // Array of pointers to hits of track

public:
   TpcLheCMTrack();                   // constructor
   TpcLheCMTrack(Int_t tracknumber);  // constructor with tracknumber
   TpcLheCMTrack(Int_t tn, Int_t nh); //

   virtual ~TpcLheCMTrack(); // destructor

   virtual void Print();

   void AddPoint(TpcLheCMPoint *point, Bool_t bward); // adds a point
   void Copy(const TpcLheCMTrack *src);               //

   // getters

   TObjArray *GetCMHits() const { return fMappingHits; }
   //  Double_t   GetChi2Bend()         const { return fChi2Bend; }
   //  Double_t   GetChi2Deep()         const { return fChi2Deep; }
   Int_t GetNumberOfPoints() const { return fNMapHits; }

   // setters
   void SetTrackID(Int_t number);
   //  void   SetChi2Bend(Double_t f)      { fChi2Bend = f; }
   //  void   SetChi2Deep(Double_t f)      { fChi2Deep = f; }
   void SetPointsUsage();
   void Clear();

   ClassDef(TpcLheCMTrack, 1) //
};

#endif
