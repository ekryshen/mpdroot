#ifndef MPDSTSGEO_H
#define MPDSTSGEO_H

#include "FairGeoSet.h"

class MpdStsGeo : public FairGeoSet {

protected:
   char modName[20]; // name of module
   char eleName[20]; // substring for elements in module

public:
   MpdStsGeo();
   ~MpdStsGeo() {}
   const char *getModuleName(Int_t);
   const char *getEleName(Int_t);

   ClassDef(MpdStsGeo, 0) // Class for Sts
};

#endif /* !MPDSTSGEO_H */
