/*
 * MpdFloorGeo.h
 *
 *  Created on: 21 maj 2019
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef Floor_Floor_MPDFloorGEO_H_
#define Floor_Floor_MPDFloorGEO_H_
#include "FairGeoSet.h"                 // for FairGeoSet

#include "Rtypes.h"                     // for Int_t, etc

#include "TString.h"
class MpdFloorGeo : public FairGeoSet{
	  protected:
    char modName[20];  // name of module
    char eleName[20]; // substring for elements in module
public:
	MpdFloorGeo();
	const char* getModuleName(Int_t );
	const char* getEleName(Int_t);
	inline Int_t getModNumInMod(const TString& name);
	virtual ~MpdFloorGeo();
	ClassDef(MpdFloorGeo,1)

};

inline Int_t MpdFloorGeo::getModNumInMod(const TString& name)
{
  /** returns the module index from module name
   ?? in name[??] has to be the length of the detector name in the
   .geo file. For example if all nodes in this file starts with
   tutdet ?? has to be 6.
  */
  return static_cast<Int_t>((name[5]-'0')-1); //
}

#endif /* Floor_Floor_MPDFloorGEO_H_ */
