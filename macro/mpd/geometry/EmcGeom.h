#ifndef EmcGeom_H
#define EmcGeom_H

#include <iostream>
#include "TObject.h"
#include "TObjArray.h"
#include "TMath.h"


#include "mpdshape.class.C"
#include "../../../detectors/emc/emc.h"
#include <vector>

using namespace std;

class EmcGeom : public TObject  
{

public:

  EmcGeom();
  EmcGeom(ofstream* f) {
    fGeoFile = f;
  }
     
  virtual ~EmcGeom();
     
  void BuildEMC();
 // int build_wall (TString w);  
 // void FieldCage(Double_t z); 
  void BuildSensVolume();
 //void BuildEC();
  private:

 ofstream* fGeoFile ;    

 ClassDef(EmcGeom,1);
};

#endif