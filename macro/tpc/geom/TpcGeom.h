#ifndef TpcGeom_H
#define TpcGeom_H

#include <iostream>
#include "TObject.h"
#include "TObjArray.h"
#include "TMath.h"

#include "tpc.h"
#include "mpdshape.class.C"


#include <vector>

using namespace std;

//class MpdShape;

class TpcGeom : public TObject  
{

public:

  TpcGeom();
  TpcGeom(ofstream* f) {
    fGeoFile = f;
  }
    virtual ~TpcGeom();
    
    /// clear lists, free pointers, etc., after read from / write to file
  //  virtual void Clear(Option_t *option ="");

  void BuildTPC();
  int build_wall (TString w);
  int build_membrane();
  void BuildSensVolume();
  void FieldCage(Double_t z);
  void BuildSensVolume_PGON();
  void padpl_layers(Double_t x_small, Double_t x_large,
		    Double_t yWidth, TList *layers);
  void BuildPadPlane();
  
  void BuildEC();
  void BuildRibs(Double_t zWidth);
  
  
private:

 ofstream* fGeoFile ;

    ClassDef(TpcGeom,1);
};

#endif
