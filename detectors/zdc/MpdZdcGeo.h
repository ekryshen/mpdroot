/*************************************************************************************
 *
 *         Class MpdZdcGeo
 *
 *  Adopted for MPD by:   Elena Litvinenko
 *  e-mail:   litvin@nf.jinr.ru
 *  Version:  8-Apr-2008
 *
 ************************************************************************************/

#ifndef MPDGEOZDC_H
#define MPDGEOZDC_H

#include "FairGeoSet.h"

class MpdZdcGeo : public FairGeoSet {
protected:
   char modName[20]; // name of module
   char eleName[20]; // substring for elements in module

public:
   MpdZdcGeo();
   ~MpdZdcGeo() {}
   const char *getModuleName(Int_t);
   const char *getEleName(Int_t);

   ClassDef(MpdZdcGeo, 0) // Class for Hyp
};
#endif /* !MPDGEOZDC_H */
