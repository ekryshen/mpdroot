//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_ETOF_GEO_H
#define __MPD_ETOF_GEO_H 1

#include "FairGeoSet.h"
//------------------------------------------------------------------------------------------------------------------------
class MpdEtofGeo : public FairGeoSet {

protected:
   char modName[20]; // name of module
   char eleName[20]; // substring for elements in module

public:
   MpdEtofGeo();
   ~MpdEtofGeo(){};

   const char *getModuleName(Int_t);
   const char *getEleName(Int_t);

   ClassDef(MpdEtofGeo, 0);
};
//------------------------------------------------------------------------------------------------------------------------
#endif
