// -------------------------------------------------------------------------
//                            MpdRegion header file                    -----
//        Created 23/07/13  by P. Batyuk (MPD) <batyuk@jinr.ru>        -----
//                          from PndRegion (PNDROOT)                   -----
// -------------------------------------------------------------------------

#ifndef MPDREGION_H
#define MPDREGION_H 1
#include "TObject.h"

class MpdRegion : public TObject {
  
 public:
  
  /** Default constructor **/
  MpdRegion( Double_t Zmin, Double_t Zmax);
  /** Destructor **/
  virtual ~MpdRegion();
  Bool_t IsInside(Double_t Z);
  ClassDef(MpdRegion,1) 
   
    protected:
  Double_t fZmin;
  Double_t fZmax;
};

#endif



