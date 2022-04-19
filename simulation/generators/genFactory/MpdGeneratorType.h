#ifndef MPDGENERATORTYPE_H
#define MPDGENERATORTYPE_H
#include <TString.h>

enum class MpdGeneratorType { URQMD, VHLLE, FLUID, PART, ION, BOX, HSD, LAQGSM, CUSTOM = 1000 };

struct MpdGenerator {
   MpdGeneratorType genType;
   TString          inFile;
   Int_t            startEvent;
   Int_t            nEvents;
};
#endif
