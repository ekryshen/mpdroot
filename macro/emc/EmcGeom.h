#ifndef EmcGeom_H
#define EmcGeom_H

#include <iostream>
#include <vector>

#include "TObject.h"
#include "TObjArray.h"
#include "TMath.h"

class EmcGeom : public TObject {

public:
   EmcGeom();
   EmcGeom(ofstream *f) { fGeoFile = f; }

   virtual ~EmcGeom();

   void BuildEMC();
   // int build_wall (TString w);
   // void FieldCage(Double_t z);
   void BuildSensVolume();
   // void BuildEC();
private:
   ofstream *fGeoFile;

   ClassDef(EmcGeom, 1);
};

#endif
