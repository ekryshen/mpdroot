#include "mpdshape.class.C"
#include "EmcGeom.h"
class Mpdshape;


  ofstream* f ;
  std::ostringstream points, position, rotation;

//____________________________________________________________________
Int_t EmcGeom() {
 
  gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C");
  mpdloadlibs(0,0); 
  
//  gROOT->ProcessLine(".L mpdshape.class.C");
//  gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/geometry/mpdshape.class.C");
  
 // gROOT->LoadClass("Mpdshape");

  const char* filename = "emc_tr_400_3.geo";

  f = new ofstream(filename, ios::out | ios::trunc);

  gROOT->ProcessLine(".L EmcGeom_lead_1_n.cxx");

  EmcGeom* emc = new EmcGeom(f);

  emc->BuildEMC();
  
  f->close(); 

}
