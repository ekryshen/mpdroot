#ifndef MPDBBCGEOPAR_H
#define MPDBBCGEOPAR_H

#include "FairParGenericSet.h"
#include "TH1F.h"

class MpdBbcGeoPar : public FairParGenericSet {
public:
   TObjArray *fGeoSensNodes; // List of FairGeoNodes for sensitive volumes
   TObjArray *fGeoPassNodes; // List of FairGeoNodes for sensitive volumes

   MpdBbcGeoPar(const char *name = "MpdBbcGeoPar", const char *title = "Bbc Geometry Parameters",
                const char *context = "TestDefaultContext");
   ~MpdBbcGeoPar(void);
   void       clear(void);
   void       putParams(FairParamList *);
   Bool_t     getParams(FairParamList *);
   TObjArray *GetGeoSensitiveNodes() { return fGeoSensNodes; }
   TObjArray *GetGeoPassiveNodes() { return fGeoPassNodes; }

   ClassDef(MpdBbcGeoPar, 1);
};

#endif /* !MPDBBCGEOPAR_H */
