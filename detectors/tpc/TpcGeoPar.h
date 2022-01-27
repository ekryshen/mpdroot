#ifndef TPCGEOPAR_H
#define TPCGEOPAR_H

#include "FairParGenericSet.h"
#include "TH1F.h"

class TpcGeoPar : public FairParGenericSet {
public:
   TObjArray *fGeoSensNodes; // List of FairGeoNodes for sensitive volumes
   TObjArray *fGeoPassNodes; // List of FairGeoNodes for sensitive volumes

   TpcGeoPar(const char *name = "TpcGeoPar", const char *title = "Tpc Geometry Parameters",
             const char *context = "TpcDefaultContext");
   ~TpcGeoPar(void);
   void       clear(void);
   void       putParams(FairParamList *);
   Bool_t     getParams(FairParamList *);
   TObjArray *GetGeoSensitiveNodes() { return fGeoSensNodes; }
   TObjArray *GetGeoPassiveNodes() { return fGeoPassNodes; }

   ClassDef(TpcGeoPar, 1);
};
#endif
