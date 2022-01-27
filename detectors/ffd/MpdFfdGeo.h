#ifndef MPDFFDGEO_H
#define MPDFFDGEO_H

#include "FairGeoSet.h"

class MpdFfdGeo : public FairGeoSet {

protected:
   char modName[20]; // name of module
   char eleName[20]; // substring for elements in module

public:
   MpdFfdGeo();
   ~MpdFfdGeo() {}
   const char *getModuleName(Int_t);
   const char *getEleName(Int_t);

   ClassDef(MpdFfdGeo, 0); // Class for Ffd
};

#endif /* !MPDFFDGEO_H */
