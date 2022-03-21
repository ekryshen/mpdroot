#ifndef CRADLE_H
#define CRADLE_H

#include "FairModule.h"
#include "Rtypes.h"

class FairCradle : public FairModule {
public:
   FairCradle(const char *name, const char *Title = "Cradle");
   FairCradle();
   virtual ~FairCradle();
   virtual void ConstructGeometry();
   ClassDef(FairCradle, 0)
};

#endif // CRADLE_H
