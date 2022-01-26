#ifndef BMDGEOPAR_H
#define BMDGEOPAR_H

#include "FairParGenericSet.h"
#include "TH1F.h"

class BmdGeoPar : public FairParGenericSet {
public:
  TObjArray            *fGeoSensNodes; // List of FairGeoNodes for sensitive volumes
  TObjArray            *fGeoPassNodes; // List of FairGeoNodes for sensitive volumes

  BmdGeoPar(const char* name="BmdGeoPar",
	    const char* title="Bmd Geometry Parameters",
	    const char* context="BmdDefaultContext");
  ~BmdGeoPar(void);
  void clear(void);
  void putParams(FairParamList*);
  Bool_t getParams(FairParamList*);
  TObjArray             *GetGeoSensitiveNodes(){return fGeoSensNodes;}
  TObjArray             *GetGeoPassiveNodes(){return fGeoPassNodes;}

  ClassDef(BmdGeoPar,1)
};

#endif
