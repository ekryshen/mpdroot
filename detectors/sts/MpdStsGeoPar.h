#ifndef MPDSTSGEOPAR_H
#define MPDSTSGEOPAR_H

#include "FairParGenericSet.h"
#include "TH1F.h"

class MpdStsGeoPar : public FairParGenericSet {
public:
   TObjArray *fGeoSensNodes; // List of FairGeoNodes for sensitive volumes
   TObjArray *fGeoPassNodes; // List of FairGeoNodes for sensitive volumes

   MpdStsGeoPar(const char *name = "MpdStsGeoPar", const char *title = "Sts Geometry Parameters",
                const char *context = "TestDefaultContext");
   ~MpdStsGeoPar(void);
   void       clear(void);
   void       putParams(FairParamList *);
   Bool_t     getParams(FairParamList *);
   TObjArray *GetGeoSensitiveNodes() { return fGeoSensNodes; }
   TObjArray *GetGeoPassiveNodes() { return fGeoPassNodes; }

   ClassDef(MpdStsGeoPar, 1);
};

#endif /* !MPDSTSGEOPAR_H */
