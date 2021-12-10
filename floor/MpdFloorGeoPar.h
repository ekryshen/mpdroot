/*
 * MpdFloorPar.h
 *
 *  Created on: 21 maj 2019
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef Floor_Floor_MPDFloorGEOPAR_H_
#define Floor_Floor_MPDFloorGEOPAR_H_

#include "FairParGenericSet.h"          // for FairParGenericSet
#include "Rtypes.h"

class TObjArray;
class FairParamList;

class MpdFloorGeoPar : public FairParGenericSet{
public:

  /** List of FairGeoNodes for sensitive  volumes */
  TObjArray*      fGeoSensNodes;

  /** List of FairGeoNodes for sensitive  volumes */
  TObjArray*      fGeoPassNodes;

  MpdFloorGeoPar(const char* name="MpdFloorGeoPar",
                         const char* title="MpdFloor Geometry Parameters",
                         const char* context="TestDefaultContext");
  ~MpdFloorGeoPar(void);
  void clear(void);
  void putParams(FairParamList*);
  Bool_t getParams(FairParamList*);
  TObjArray* GetGeoSensitiveNodes() {return fGeoSensNodes;}
  TObjArray* GetGeoPassiveNodes()   {return fGeoPassNodes;}

private:
  MpdFloorGeoPar(const MpdFloorGeoPar&);
  MpdFloorGeoPar& operator=(const MpdFloorGeoPar&);

ClassDef(MpdFloorGeoPar,1)
};

#endif /* Floor_Floor_MPDFloorGEOPAR_H_ */
