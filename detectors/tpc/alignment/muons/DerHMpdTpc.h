#ifndef DERHMPDTPC_H
#define DERHMPDTPC_H

/// \ingroup
/// \class DerHMpdTpc
/// \brief Partial derivatives for TPC alignment
///
/// \author Valentin Kuzmin, SINP of Moscow State University

using namespace std;
#include <iostream>
#include <TObject.h>
#include <TVector3.h>
#include <TRotation.h>
#include <TMath.h>
#include <miniTpcAlignment.h>

// Collaborating Class Headers -------

#include "TpcSectorGeoAlignmentVAK.h"

// Collaborating Class Declarations --

class DerHMpdTpc : public TObject
/*
  The MpdTpcMiniMC class calculates partial derivatives for
  Millipede allignment procedure.
  It is the line track model with 5 parameters

*/
{

public:
   ///< Constructor
   DerHMpdTpc(BaseTpcSectorGeo &tpcGeo)
   {
      fTpcSecGeo = dynamic_cast<TpcSectorGeoAlignmentVAK *>(&tpcGeo);
      if (!fTpcSecGeo)
         Fatal("DerHMpdTpc::DerHMpdTpc", " !!! Wrong geometry type !!! "); // pointer to TPC geometry
                                                                           // to get TPC Sector alignment
      Init();
   }

   ///< Destructor
   virtual ~DerHMpdTpc() { ; } ///< Destructor

   void Init();
   void SetMagFld(TRotation g2b, TRotation b2g)
   {
      G2B = g2b;
      B2G = b2g;
   };                                           // GCS<->MFCS
   void SetLocPar(Double_t *qt);                // set local parameters
   void SetHitTrack(Int_t isec, Double_t *hit); // set hits in LCS
   void GetLocDer(Double_t *dFdq);              // get local derivatives
   void GetGloDer(Int_t isec, Double_t *dFdp);  // get global derivatives
   //  void GetGloDerT(Int_t isec,Double_t* dFdp);           //get global derivatives
   void     GetGlodHdp(Int_t isec, Double_t dhdp[3][6]);
   void     Getd2Fdpdp(Int_t isec, Double_t d2Fdpdp[6][6]);
   void     GetLocDer2(Double_t d2Fdqdq[6][6]);
   void     GetLocGloDer2(Int_t isec, Double_t d2Fdqdp[6][6]);
   Double_t gett0() { return t0; };

   void PrintdR(const char *title, Int_t s, Double_t R[3][3]);
   void PrintM66(const char *title, Double_t m[6][6]);
   //______________________________________________________________________

private:
   // Sector geometry
   TpcSectorGeoAlignmentVAK *fTpcSecGeo;
   // Alignment parameters (global parameters)
   TVector3  R0_gsc[24], G0_gsc[24], R0_tlsc[24], G0_tlsc[24];
   TRotation G2B, B2G;
   TRotation glo2loc[24];
   TRotation Tlo2glo[24]; // ideal sector TLCS -> GCS transformation (Z-rotation)
   TRotation Tglo2lo[24]; // GCS -> ideal sector TLCS transformation (Z-rotation)
                          // Track model parameters (local parameters)
   Double_t q[6];
   Double_t dAdp[3][3][3][24], dRdp[3][3][3][24];
   Double_t d2Adp2[3][3][3][3][24], d2Rdp2[3][3][3][3][24];
   Double_t dHdp[3][6];
   Double_t dTdq[3][6];
   Double_t d2Tdqdq[3][6][6];
   // Unit vector of the track line
   Double_t Rh, q_ch, e[3], sq4, cq4, sq5, cq5;
   // Hit & Track coordinates
   Double_t t0, h[3], ht[3], hs[3], T[3], Tt[3];

   // -----  ----- methods --------------

   ClassDef(DerHMpdTpc, 1);
};

#endif
