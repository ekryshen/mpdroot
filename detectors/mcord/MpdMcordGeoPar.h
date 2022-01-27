/*
 * MpdMcordPar.h
 *
 *  Created on: 21 maj 2019
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MCORD_MCORD_MPDMCORDGEOPAR_H_
#define MCORD_MCORD_MPDMCORDGEOPAR_H_
#include "FairParGenericSet.h" // for FairParGenericSet
#include "Rtypes.h"

class TObjArray;
class FairParamList;

class MpdMcordGeoPar : public FairParGenericSet {
public:
   /** List of FairGeoNodes for sensitive  volumes */
   TObjArray *fGeoSensNodes;

   /** List of FairGeoNodes for sensitive  volumes */
   TObjArray *fGeoPassNodes;

   MpdMcordGeoPar(const char *name = "MpdMcordGeoPar", const char *title = "MpdMcord Geometry Parameters",
                  const char *context = "TestDefaultContext");
   ~MpdMcordGeoPar(void);
   void       clear(void);
   void       putParams(FairParamList *);
   Bool_t     getParams(FairParamList *);
   TObjArray *GetGeoSensitiveNodes() { return fGeoSensNodes; }
   TObjArray *GetGeoPassiveNodes() { return fGeoPassNodes; }

private:
   MpdMcordGeoPar(const MpdMcordGeoPar &);
   MpdMcordGeoPar &operator=(const MpdMcordGeoPar &);

   ClassDef(MpdMcordGeoPar, 1);
};
#endif /* MCORD_MCORD_MPDMCORDGEOPAR_H_ */
