#ifndef MPDEMCDETGEO_H
#define MPDEMCDETGEO_H

#include "FairGeoSet.h"

class  MpdEmcGeo : public FairGeoSet {

 protected:
  char modName[20];  // name of module
  char eleName[20];  // substring for elements in module
 public:
  MpdEmcGeo();
  ~MpdEmcGeo() {}
  const char* getModuleName(Int_t);
  const char* getEleName(Int_t);
  inline Int_t getModNumInMod(const TString&);
  ClassDef(MpdEmcGeo,1)
};

inline Int_t MpdEmcGeo::getModNumInMod(const TString& name) {
  /** returns the module index from module name 
   ?? in name[??] has to be the length of the detector name in the
   .geo file. For example if all nodes in this file starts with
   tutdet ?? has to be 6.
  */
  return (Int_t)(name[6]-'0')-1; // 
}

#endif  /*  MPDEMCDETGEO_H */
