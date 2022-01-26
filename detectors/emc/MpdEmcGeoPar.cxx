#include "MpdEmcGeoPar.h"
#include "FairParamList.h"
#include "TObjArray.h"
#include "emc.h"
#include <iostream>

ClassImp(MpdEmcGeoPar)

MpdEmcGeoPar ::MpdEmcGeoPar(const char* name, const char* title, const char* context)
  : FairParGenericSet(name,title,context) {

  fGeoSensNodes = new TObjArray();
  fGeoPassNodes = new TObjArray();

}

MpdEmcGeoPar::MpdEmcGeoPar() {
    //TODO!!! correct parameters;
    length = emc::emc1Chamber_z * 0.1; //0.1: mm->cm
    rMax = emc::outr;
    rMin = emc::inr;
    nSec = emc::NSector;
    nRows = emc::NmodZ;
    nSupMod = emc::nSuperModInRow;
    nMod = emc::nModInSuperModByPhi * emc::nModInSuperModByZ;
    lengthOfModuleByZ = emc::base * 0.1; // 0.1: mm->cm
    lengthOfSuperModuleByZ = emc::emc1_box1_z * 0.1; // 0.1: mm->cm
    nModInSuperModByZ = emc::nModInSuperModByZ; // 0.1: mm->cm
    nModInSuperModByPhi = emc::nModInSuperModByPhi;
    angleOfSuperModule = emc::angleSuperModule;
    angleOfSector = emc::angleSector;
    angleOfModule = emc::angleModule;
}

MpdEmcGeoPar::~MpdEmcGeoPar(void) {
}

void MpdEmcGeoPar::clear(void) {
  if(fGeoSensNodes) delete fGeoSensNodes;
  if(fGeoPassNodes) delete fGeoPassNodes;
}

void MpdEmcGeoPar::putParams(FairParamList* l) {
  if (!l) return;
  l->addObject("FairGeoNodes Sensitive List", fGeoSensNodes);
  l->addObject("FairGeoNodes Passive List", fGeoPassNodes);
}

Bool_t MpdEmcGeoPar::getParams(FairParamList* l) {
  if (!l) return kFALSE;
  if (!l->fillObject("FairGeoNodes Sensitive List", fGeoSensNodes)) return kFALSE;
  if (!l->fillObject("FairGeoNodes Passive List", fGeoPassNodes)) return kFALSE;
  return kTRUE;
}