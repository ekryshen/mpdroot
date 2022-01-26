#ifndef MPDEMCDETGEOPAR_H
#define MPDEMCDETGEOPAR_H

#include "FairParGenericSet.h"

class TObjArray;
class FairParamList;

class MpdEmcGeoPar       : public FairParGenericSet {
public:
 
  /** List of FairGeoNodes for sensitive  volumes */
  TObjArray      *fGeoSensNodes; 

  /** List of FairGeoNodes for sensitive  volumes */
  TObjArray      *fGeoPassNodes; 

  MpdEmcGeoPar(const char* name, const char* title, const char* context);
  MpdEmcGeoPar();
  ~MpdEmcGeoPar(void);
  void clear(void);
  void putParams(FairParamList*);
  Bool_t getParams(FairParamList*);
  TObjArray* GetGeoSensitiveNodes() {return fGeoSensNodes;}
  TObjArray* GetGeoPassiveNodes()   {return fGeoPassNodes;}
  
  Float_t GetLength() const {return length;}
  Float_t GetRmin() const {return rMin;}
  Float_t GetRmax() const {return rMax;}
  
  UInt_t GetNsec() const {return nSec;}
  UInt_t GetNrows() const {return nRows;}
  UShort_t GetNmod() const {return nMod;}
  UInt_t GetNsupMod() const {return nSupMod;}
  
  Float_t GetLengthOfModuleByZ() const {return lengthOfModuleByZ;}
  Float_t GetAngleOfModule() const {return angleOfModule;}
  Float_t GetAngleOfSuperModule() const {return angleOfSuperModule;}
  Float_t GetAngleOfSector() const {return angleOfSector;}
  Float_t GetLengthOfSuperModuleByZ() const {return lengthOfSuperModuleByZ;}
  UInt_t GetNModInSuperModByZ() const {return nModInSuperModByZ;}
  UInt_t GetNModInSuperModByPhi() const {return nModInSuperModByPhi;}
  
private:
    
  Float_t length;
  Float_t rMin;
  Float_t rMax;
  Float_t lengthOfModuleByZ;
  Float_t lengthOfSuperModuleByZ;
  Float_t angleOfSuperModule;
  Float_t angleOfModule;
  Float_t angleOfSector;
  
  UInt_t nModInSuperModByZ;
  UInt_t nModInSuperModByPhi;
  UInt_t nSec;    // number of sectors in whole EMC (positive and negative Z)
  UInt_t nRows;   // number of rows (slices) in one sector
  UShort_t nMod;  // number of modules in one super-module (N x M)
  UInt_t nSupMod; // number of super-modules in one row (slice)

  ClassDef(MpdEmcGeoPar, 1)
};

#endif /* MPDEMCDETGEOPAR_H */
