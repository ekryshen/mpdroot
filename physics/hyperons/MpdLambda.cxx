/// \ingroup physics
/// \class MpdLambda
/// \brief Lambda reconstructed object
///
/// \author Alexander Zinchenko, LHEP JINR Dubna - 22-03-2023

#include "MpdLambda.h"

ClassImp(MpdLambda);

//-----------------------------------------------------------------------------------------------

MpdLambda::MpdLambda(Float_t imassh, Float_t ipth, Float_t ipthgen, Float_t iph, Float_t ietah, Float_t iyh,
                     Float_t ichi2h, Float_t idisth, Float_t ipath, Float_t iangle, Float_t *ietas, Float_t *ips,
                     Float_t *ipts, Float_t *ichi2s, Float_t *idcas, Float_t *ic2s, Int_t *iorigs, Int_t *iqs,
                     Int_t *ilayMx, Int_t impdg)
   : massh(imassh), pth(ipth), pthgen(ipthgen), ph(iph), etah(ietah), yh(iyh), chi2h(ichi2h), disth(idisth),
     path(ipath), angle(iangle), mpdg(impdg)
{

   for (Int_t j = 0; j < 2; ++j) {
      etas[j]  = ietas[j];
      ps[j]    = ips[j];
      pts[j]   = ipts[j];
      chi2s[j] = ichi2s[j];
      dcas[j]  = idcas[j];
      c2s[j]   = ic2s[j];
      origs[j] = iorigs[j];
      qs[j]    = iqs[j];
      layMx[j] = ilayMx[j];
   }
}

//-----------------------------------------------------------------------------------------------
