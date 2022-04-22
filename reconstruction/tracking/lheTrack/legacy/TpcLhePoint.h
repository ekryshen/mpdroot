#ifndef LHE_POINT_H
#define LHE_POINT_H

// *************************************************************************
//  Author: Oleg Rogachevsky             e-mail: rogach@sunhe.jinr.ru
//
//  store information for point
//
// Created: 1-07-07
// Modified:
//
// *************************************************************************
#include "TObject.h"
#include "TVector3.h"

class TpcLhePoint : public TObject {

private:
   TVector3 fCoord; // coordinates of point
   TVector3 fError; // errors of the coordinates

public:
   TpcLhePoint(); //

   // constructor from array of doubles
   TpcLhePoint(Double_t pos[3], Double_t err[3]);

   // constructor from doubles
   TpcLhePoint(Double_t x, Double_t y, Double_t z);
   virtual ~TpcLhePoint();

   // getters
   Double_t GetX() const { return fCoord.X(); }
   Double_t GetY() const { return fCoord.Y(); }
   Double_t GetZ() const { return fCoord.Z(); }
   Double_t GetXerr() const { return fError.X(); }
   Double_t GetYerr() const { return fError.Y(); }
   Double_t GetZerr() const { return fError.Z(); }
   TVector3 GetCoord() const { return fCoord; }
   TVector3 GetError() const { return fError; }

   // setters
   void SetX(Double_t f) { fCoord.SetX(f); }
   void SetY(Double_t f) { fCoord.SetY(f); }
   void SetZ(Double_t f) { fCoord.SetZ(f); }
   void SetXerr(Double_t f) { fError.SetX(f); }
   void SetYerr(Double_t f) { fError.SetY(f); }
   void SetZerr(Double_t f) { fError.SetZ(f); }

   ClassDef(TpcLhePoint, 1); //
};

#endif
