#ifndef MPDBBCGEO_H
#define MPDBBCGEO_H

#include "FairGeoSet.h"

class MpdBbcGeo : public FairGeoSet {

protected:
   char modName[20]; // name of module
   char eleName[20]; // substring for elements in module

public:
   MpdBbcGeo();
   ~MpdBbcGeo() {}
   const char *getModuleName(Int_t);
   const char *getEleName(Int_t);

   ClassDef(MpdBbcGeo, 0) // Class for Bbc
};

#endif /* !MPDBBCGEO_H */
