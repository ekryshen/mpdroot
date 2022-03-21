#ifndef GEOCRADLE_H
#define GEOCRADLE_H

#include "FairGeoSet.h"
#include "Rtypes.h"

class  FairGeoCradle : public FairGeoSet {
protected:
    char modName[2];  // name of module
    char eleName[2];  // substring for elements in module

public:
    FairGeoCradle();
    ~FairGeoCradle() {}
    const char* getModuleName(Int_t) {return modName;}
    const char* getEleName(Int_t) {return eleName;}
    ClassDef(FairGeoCradle,0) // Class for geometry of Cradle
};

#endif  /* !GEOCRADLE_H */
