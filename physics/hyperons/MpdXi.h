#ifndef MPDXI_H
#define MPDXI_H

/// \ingroup physics
/// \class MpdXi
/// \brief Xi reconstructed object
///
/// \author Alexander Zinchenko, LHEP JINR Dubna - 25-03-2023

#include "TObject.h"

class MpdXi {

public:
   MpdXi() { ; }
   MpdXi(Float_t imassh, Float_t ipth, Float_t ipthgen, Float_t iph, Float_t ietah, Float_t iyh, Float_t ichi2h,
         Float_t idisth, Float_t ipath, Float_t iangle, Float_t ic2pv, Float_t idca, Float_t *ietas, Float_t *ips,
         Float_t *ipts, Float_t *ichi2s, Float_t *idcas, Float_t *ic2s, Float_t imassL, Float_t ichi2L, Float_t idistL,
         Float_t ipathL, Float_t iangL, Float_t *ichi2sL, Float_t imassL1, Int_t *iorigs, Int_t *iqs, Int_t *ilayMx,
         Int_t impdg);
   virtual ~MpdXi() {}

public:
   Float_t massh, pth, pthgen, ph, etah, yh, chi2h, disth, path, angle, c2pv, dca, etas[2];
   Float_t ps[2], pts[2], chi2s[2], dcas[2], massL, chi2L, distL, pathL, angL;
   Float_t chi2sL[2], massL1, c2s[3];
   Int_t   origs[2], qs[2], layMx[2], mpdg;

   // private:

   ClassDef(MpdXi, 2);
};
#endif
