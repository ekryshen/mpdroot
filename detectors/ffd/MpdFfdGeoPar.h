#ifndef MPDFFDGEOPAR_H
#define MPDFFDGEOPAR_H

#include "FairParGenericSet.h"
#include "TH1F.h"

class MpdFfdGeoPar : public FairParGenericSet {
public:
   TObjArray *fGeoSensNodes; // List of FairGeoNodes for sensitive volumes
   TObjArray *fGeoPassNodes; // List of FairGeoNodes for sensitive volumes

   MpdFfdGeoPar(const char *name = "MpdFfdGeoPar", const char *title = "Ffd Geometry Parameters",
                const char *context = "TestDefaultContext");
   ~MpdFfdGeoPar(void);
   void       clear(void);
   void       putParams(FairParamList *);
   Bool_t     getParams(FairParamList *);
   TObjArray *GetGeoSensitiveNodes() { return fGeoSensNodes; }
   TObjArray *GetGeoPassiveNodes() { return fGeoPassNodes; }

   ClassDef(MpdFfdGeoPar, 1);
};

#endif /* !MPDFFDGEOPAR_H */
