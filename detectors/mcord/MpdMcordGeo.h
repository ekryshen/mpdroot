/*
 * MpdMcordGeo.h
 *
 *  Created on: 21 maj 2019
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MCORD_MCORD_MPDMCORDGEO_H_
#define MCORD_MCORD_MPDMCORDGEO_H_
#include "FairGeoSet.h" // for FairGeoSet

#include "Rtypes.h" // for Int_t, etc

#include "TString.h"
class MpdMcordGeo : public FairGeoSet {
protected:
   char modName[20]; // name of module
   char eleName[20]; // substring for elements in module
public:
   MpdMcordGeo();
   const char * getModuleName(Int_t);
   const char * getEleName(Int_t);
   inline Int_t getModNumInMod(const TString &name);
   virtual ~MpdMcordGeo();
   ClassDef(MpdMcordGeo, 1);
};

inline Int_t MpdMcordGeo::getModNumInMod(const TString &name)
{
   /** returns the module index from module name
    ?? in name[??] has to be the length of the detector name in the
    .geo file. For example if all nodes in this file starts with
    tutdet ?? has to be 6.
   */
   return static_cast<Int_t>((name[5] - '0') - 1); //
}

#endif /* MCORD_MCORD_MPDMCORDGEO_H_ */
