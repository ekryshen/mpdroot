/*
 * File:   MpdV0.h
 * Author: niki
 *
 * Created on September 6, 2013, 7:10 PM
 */
#ifndef MPD_V0_H
#define MPD_V0_H

//#include "MpdVertex.h"
#include "TROOT.h"
#include "TNamed.h"
#include "Rtypes.h"
#include "TVector3.h"
#include "TClass.h"
#include "TMath.h"

class MpdV0 : public TNamed {
public:
   MpdV0();
   MpdV0(const char *name, const char *title);
   MpdV0(const char *name, const char *title, Double_t vtx[3], Double_t chi2, Double_t p[3], Double_t DCA,
         Double_t cosa, Double_t mass, Double_t r[3], Int_t B, Int_t id1, Int_t id2, Int_t id0);
   MpdV0 &operator=(const MpdV0 &V0);

   virtual ~MpdV0();

   void Print();
   void SetV0(Double_t vtx[3], Double_t chi2, Double_t p[3], Double_t DCA, Double_t cosa, Double_t mass, Double_t r[3],
              Int_t B, Int_t id1, Int_t id2, Int_t id0);

   // void SetP(TVector3 P){fp = P;}
   Double_t *GetPxPyPz(Double_t p[3]);
   Double_t *GetXYZ(Double_t vtx[3]);
   Double_t *GetRxRyRz(Double_t r[3]);
   Double_t  GetP();
   Double_t  GetDist();
   Double_t  GetMass() { return fmass; }
   Double_t  GetBarionCharge() { return fB; }
   Int_t     GetID0() { return fID0; }
   Int_t     GetID1() { return fID1; }
   Int_t     GetID2() { return fID2; }
   Double_t  GetCosPA() { return fcosa; }
   Double_t  GetDCA() { return fDCA; }
   Double_t  GetChi2() { return fchi2; }

private:
   Int_t    fID0;
   Int_t    fID1;  // mc ID of first track
   Int_t    fID2;  //   of second
   Int_t    fB;    // barion charge
   Double_t fchi2; // chi square of vertex fit
   Double_t fmass; // mass; yup...
   Double_t fcosa; // pointing angle between R-v0 and P
   Double_t fDCA;  // daughters' Distance of Closest Approach
                   // Momentum
   Double_t fPx;
   Double_t fPy;
   Double_t fPz;
   Double_t fRx;
   Double_t fRy;
   Double_t fRz; // R=pvert-vtx or flight - distance between primary and v0
   Double_t fX;
   Double_t fY;
   Double_t fZ; // vtx(x,y,z)

   // const TMatrixD &covMat;       //Covariance Matrix (symmetric; 3x3)
   ClassDef(MpdV0, 2);
};
#endif /* MPDV0_H */
