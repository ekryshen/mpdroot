//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_TOF_GEO_H
#define __MPD_TOF_GEO_H

#include "FairGeoSet.h"
//------------------------------------------------------------------------------------------------------------------------
class MpdTofGeo : public FairGeoSet {
protected:
   char modName[20]; // name of module
   char eleName[20]; // substring for elements in module

public:
   MpdTofGeo();
   ~MpdTofGeo(){};

   const char * getModuleName(Int_t);
   const char * getEleName(Int_t);
   inline Int_t getModNumInMod(const TString &);

   ClassDef(MpdTofGeo, 0); // Class for Tof
};

//------------------------------------------------------------------------------------------------------------------------
inline Int_t MpdTofGeo::getModNumInMod(const TString &name)
{
   // returns the module index from module name
   return (Int_t)(name[3] - '0') - 1;
}
//------------------------------------------------------------------------------------------------------------------------
#endif // #ifndef __MPD_TOF_GEO_H
