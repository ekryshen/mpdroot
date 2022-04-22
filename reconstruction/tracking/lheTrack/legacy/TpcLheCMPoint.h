#ifndef LHE_CM_POINT_H
#define LHE_CM_POINT_H

// *************************************************************************
//  Author: Oleg Rogachevsky             e-mail: rogach@sunhe.jinr.ru
//
// conformal mapping point class
//
// Created: 1-07-07
// Modified:
//
// *************************************************************************
#include "TpcLheHit.h"
#include "TpcLhePoint.h"

class TpcLheCMPoint : public TpcLheHit {

private:
   // Interaction point
   Double_t fXt; // x-value of the interaction point
   Double_t fYt; // y-value of the interaction point
   Double_t fZt; // z-value of the interaction point

   Double_t fXterr; // error of fXt
   Double_t fYterr; // error of fYt
   Double_t fZterr; // error of fZt

   // conformal mapping coordinates
   Double_t fXprime; // transformed x
   Double_t fYprime; // transformed y

   Double_t fXprimeerr; // error of fXprime
   Double_t fYprimeerr; // error of fYprime

   // coordinates with respect to the vertex

   // cartesian coordinates
   Double_t fXv; // x with respect to vertex
   Double_t fYv; // y with respect to vertex
   Double_t fZv; // z with respect to vertex

   Double_t fXverr; // error of fXv
   Double_t fYverr; // error of fYv
   Double_t fZverr; // error of fZv

   // spherical coordinates
   Double_t fPhi;   // angle phi
   Double_t fTheta; // angle eta

   // distances to the conformal mapping fits
   Double_t fCircleDist; // distance from the fitted circle
   Double_t fDeepDist;   // distance from the line fit in YZ plane

public:
   TpcLheCMPoint();                 //
   TpcLheCMPoint(TpcLheHit *point); // from hit

   virtual ~TpcLheCMPoint(); //

   // getters
   Double_t GetCircleDist() const { return fCircleDist; }
   Double_t GetDeepDist() const { return fDeepDist; }

   Double_t GetPhi() const { return fPhi; }
   Double_t GetTheta() const { return fTheta; }

   Double_t GetXprime() const { return fXprime; }
   Double_t GetYprime() const { return fYprime; }
   Double_t GetXprimeerr() const { return fXprimeerr; }
   Double_t GetYprimeerr() const { return fYprimeerr; }

   Double_t GetXt() const { return fXt; }
   Double_t GetYt() const { return fYt; }
   Double_t GetZt() const { return fZt; }
   Double_t GetXterr() const { return fXterr; }
   Double_t GetYterr() const { return fYterr; }
   Double_t GetZterr() const { return fZterr; }

   Double_t GetXv() const { return fXv; }
   Double_t GetYv() const { return fYv; }
   Double_t GetZv() const { return fZv; }
   Double_t GetXverr() const { return fXverr; }
   Double_t GetYverr() const { return fYverr; }
   Double_t GetZverr() const { return fZverr; }

   // setters
   void SetDist(Double_t c, Double_t l)
   {
      fCircleDist = c;
      fDeepDist   = l;
   }
   void SetCircleDist(Double_t f) { fCircleDist = f; }
   void SetDeepDist(Double_t f) { fDeepDist = f; }

   void SetPhi(Double_t f) { fPhi = f; }
   void SetTheta(Double_t f) { fTheta = f; }

   void SetXt(Double_t f) { fXt = f; }
   void SetYt(Double_t f) { fYt = f; }
   void SetZt(Double_t f) { fZt = f; }
   void SetXterr(Double_t f) { fXterr = f; }
   void SetYterr(Double_t f) { fYterr = f; }
   void SetZterr(Double_t f) { fZterr = f; }

   void SetXv(Double_t f) { fXv = f; }
   void SetYv(Double_t f) { fYv = f; }
   void SetZv(Double_t f) { fZv = f; }
   void SetXverr(Double_t f) { fXverr = f; }
   void SetYverr(Double_t f) { fYverr = f; }
   void SetZverr(Double_t f) { fZverr = f; }

   void Setup(TpcLhePoint *vertex); // does the usual setup in the right order
   void SetAngles();                // calculate spherical angles and set values
   void Print();
   void SetIntPoint(const Double_t in_x = 0., const Double_t in_y = 0., const Double_t in_z = 0.,
                    const Double_t in_x_err = 0., const Double_t in_y_err = 0., const Double_t in_z_err = 0.);
   void SetShiftedCoord();                     // set shifted coordinates
   void SetConfCoord();                        // conformal mapping of coordinates
   void SetAllCoord(const TpcLheCMPoint *hit); //

   ClassDef(TpcLheCMPoint, 1); //
};

#endif
