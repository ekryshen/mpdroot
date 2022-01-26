//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_ETOF_GEO_PAR_H
#define __MPD_ETOF_GEO_PAR_H 1

#include "FairParGenericSet.h"
#include "TH1F.h"
#include "TObjArray.h"
//------------------------------------------------------------------------------------------------------------------------
class MpdEtofGeoPar : public FairParGenericSet {
   TObjArray *fGeoSensNodes; // List of FairGeoNodes for sensitive volumes
   TObjArray *fGeoPassNodes; // List of FairGeoNodes for passive volumes

public:
   MpdEtofGeoPar(const char *name = "MpdEtofGeoPar", const char *title = "ETof Geometry Parameters",
                 const char *context = "TestDefaultContext");
   ~MpdEtofGeoPar(void);

   void       clear(void);
   void       putParams(FairParamList *);
   Bool_t     getParams(FairParamList *);
   TObjArray *GetGeoSensitiveNodes() { return fGeoSensNodes; }
   TObjArray *GetGeoPassiveNodes() { return fGeoPassNodes; }

   ClassDef(MpdEtofGeoPar, 1)
};
//------------------------------------------------------------------------------------------------------------------------
#endif
