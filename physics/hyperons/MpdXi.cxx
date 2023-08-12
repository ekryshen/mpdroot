/// \ingroup physics
/// \class MpdXi
/// \brief Xi reconstructed object
///
/// \author Alexander Zinchenko, LHEP JINR Dubna - 25-03-2023

#include "MpdXi.h"

ClassImp(MpdXi);

//-----------------------------------------------------------------------------------------------

MpdXi::MpdXi(Float_t imassh, Float_t ipth, Float_t iph, Float_t ietah, Float_t iyh, Float_t ichi2h, Float_t idisth,
             Float_t ipath, Float_t iangle, Float_t ic2pv, Float_t idca, Float_t *ietas, Float_t *ips, Float_t *ipts,
             Float_t *ichi2s, Float_t *idcas, Float_t *ic2s, Float_t imassL, Float_t ichi2L, Float_t idistL,
             Float_t ipathL, Float_t iangL, Float_t *ichi2sL, Float_t imassL1, Int_t *iorigs, Int_t *iqs, Int_t *ilayMx,
             Int_t impdg)
   : massh(imassh), pth(ipth), ph(iph), etah(ietah), yh(iyh), chi2h(ichi2h), disth(idisth), path(ipath), angle(iangle),
     c2pv(ic2pv), dca(idca), massL(imassL), chi2L(ichi2L), distL(idistL), pathL(ipathL), angL(iangL), massL1(imassL1),
     mpdg(impdg)
{

   for (Int_t j = 0; j < 2; ++j) {
      etas[j]   = ietas[j];
      ps[j]     = ips[j];
      pts[j]    = ipts[j];
      chi2s[j]  = ichi2s[j];
      chi2sL[j] = ichi2sL[j];
      dcas[j]   = idcas[j];
      origs[j]  = iorigs[j];
      qs[j]     = iqs[j];
      layMx[j]  = ilayMx[j];
   }
   for (int j = 0; j < 3; ++j) c2s[j] = ic2s[j];
}

//-----------------------------------------------------------------------------------------------
