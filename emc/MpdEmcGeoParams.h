#ifndef MPDEMCDETGEOPARAMS_H
#define MPDEMCDETGEOPARAMS_H

#include "FairParGenericSet.h"
#include "TMath.h"
#include <map>

class TObjArray;
class FairParamList;

using namespace std;
using namespace TMath;

class MpdEmcGeoParams       : public FairParGenericSet {

public:
 
  /** List of FairGeoNodes for sensitive  volumes */
  TObjArray      *fGeoSensNodes; 
  /** List of FairGeoNodes for sensitive  volumes */
  TObjArray      *fGeoPassNodes; 

  MpdEmcGeoParams(const char* name, const char* title, const char* context);
  MpdEmcGeoParams();
  ~MpdEmcGeoParams(void);

  static MpdEmcGeoParams* GetInstance() {
    if (!fGeomParams) fGeomParams = new MpdEmcGeoParams();
    return fGeomParams;
  }


  void clear(void);
  void putParams(FairParamList*);
  Bool_t getParams(FairParamList*);
  TObjArray* GetGeoSensitiveNodes() {return fGeoSensNodes;}
  TObjArray* GetGeoPassiveNodes()   {return fGeoPassNodes;}
  
  const Double_t GetLength() const {return length;}
  const Double_t GetRmin() const {return rMin;}
  const Double_t GetRmax() const {return rMax;}  
  const Double_t GetXTower(Int_t iID){return xTower[pointID[iID]];}
  const Double_t GetYTower(Int_t iID){return yTower[pointID[iID]];}
  const Double_t GetZTower(Int_t iID){return zTower[pointID[iID]];}
  const Double_t GetPhiTower(Int_t iID){return phiTower[pointID[iID]];} 
  const Double_t GetThetaTower(Int_t iID){return thetaTower[pointID[iID]];} 
 
  vector<Double_t> GetPhiRow(){return phiRow;}
  vector<Double_t> GetZCenterBox() {return zTower;}
  vector<Double_t> GetThetaBox() {return thetaTower;}
  vector<Double_t> GetPhiBox() {return phiTower;}
  vector<Double_t> GetRhoCenterBox() {return rhoTower;}

  const UInt_t GetNsec() const {return fSectorNumber;}
  const UInt_t GetNrows() const {return fTowerNumberXY;}
  const UInt_t GetNmod() const {return fTowerNumberZ;}

  const Int_t GetCrateNumber(){return fCrateNumber;}
  const Int_t GetModuleTypes(){return fModuleTypes;}
  const Int_t GetTowerNumber(){return fTowerNumber;}
  const Double_t GetLengthBox(){return fTowerLength;}

private:

  static MpdEmcGeoParams* fGeomParams;

  Double_t length; // Total EMC length
  Double_t rMin; // Barell minimal radius 
  Double_t rMax; // Barell maximal radius

  Int_t fSectorNumber; // number of sectors
  Int_t fCrateNumber; // number of module lines in sector
  Int_t fModuleTypes; // number of module types 
  Int_t fTowerNumber; // total number of towers in ECal
  Int_t fTowerNumberXY; // number of towers in module (XY-plane)
  Int_t fTowerNumberZ; // number of towers in module along Z-axis
  Double_t fTowerLength; // tower length  

  map<Int_t, Int_t> pointID;
  vector<Double_t> xTower, yTower, zTower, rhoTower, phiRow; 
  vector<Double_t> thetaTower, phiTower; 


  ClassDef(MpdEmcGeoParams, 1)
};

#endif /* MPDEMCDETGEOPARAMS_H */
