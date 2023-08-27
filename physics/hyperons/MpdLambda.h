#ifndef MPDLAMBDA_H
#define MPDLAMBDA_H

/// \ingroup physics
/// \class MpdLambda
/// \brief Lambda reconstructed object
///
/// \author Alexander Zinchenko, LHEP JINR Dubna - 22-03-2023

#include "TObject.h"

class MpdLambda {
   // class MpdLambda : public TObject {

public:
   MpdLambda() { ; }
   MpdLambda(Float_t imassh, Float_t ipth, Float_t ipthgen, Float_t iph, Float_t ietah, Float_t iyh, Float_t ichi2h,
             Float_t idisth, Float_t ipath, Float_t iangle, Float_t *ietas, Float_t *ips, Float_t *ipts,
             Float_t *ichi2s, Float_t *idcas, Float_t *ic2s, Int_t *iorigs, Int_t *iqs, Int_t *ilayMx, Int_t impdg);
   virtual ~MpdLambda() {}

public:
   Float_t massh, pth, pthgen, ph, etah, yh, chi2h, disth, path, angle, etas[2], ps[2], pts[2];
   Float_t chi2s[2], c2s[2], dcas[2];
   Int_t   origs[2], qs[2], layMx[2], mpdg;

   // private:

   ClassDef(MpdLambda, 2);
};
#endif
