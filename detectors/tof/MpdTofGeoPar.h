//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_TOF_GEO_PAR_H
#define __MPD_TOF_GEO_PAR_H

#include "FairParGenericSet.h"
#include "TH1F.h"
//------------------------------------------------------------------------------------------------------------------------
class MpdTofGeoPar : public FairParGenericSet {
public:
   TObjArray *fGeoSensNodes; // List of FairGeoNodes for sensitive volumes
   TObjArray *fGeoPassNodes; // List of FairGeoNodes for sensitive volumes

   MpdTofGeoPar(const char *name = "MpdTofGeoPar", const char *title = "Tof Geometry Parameters",
                const char *context = "TestDefaultContext");
   ~MpdTofGeoPar(void);

   void       clear(void);
   void       putParams(FairParamList *);
   Bool_t     getParams(FairParamList *);
   TObjArray *GetGeoSensitiveNodes() { return fGeoSensNodes; }
   TObjArray *GetGeoPassiveNodes() { return fGeoPassNodes; }

   ClassDef(MpdTofGeoPar, 1)
};

//------------------------------------------------------------------------------------------------------------------------
#endif // #ifndef __MPD_TOF_GEO_PAR_H
