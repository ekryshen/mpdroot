#include "FairCradle.h"

FairCradle::FairCradle() : FairModule("FairCradle", "") {}

FairCradle::FairCradle(const char *name, const char *title) : FairModule(name, title) {}

FairCradle::~FairCradle() {}

void FairCradle::ConstructGeometry()
{
   TString fileName = GetGeometryFileName();
   if (fileName.EndsWith(".geo")) {
      ConstructASCIIGeometry();
   } else if (fileName.EndsWith(".root")) {
      ConstructRootGeometry();
   } else {
      std::cout << "Geometry format not supported " << std::endl;
   }
}

ClassImp(FairCradle);
