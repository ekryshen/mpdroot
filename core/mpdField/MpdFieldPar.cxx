// -------------------------------------------------------------------------
//                            MpdFieldPar source file                  -----
//        Created 23/07/13  by P. Batyuk (MPD) <batyuk@jinr.ru>        -----
//                          from MpdFieldPar (PNDROOT)                 -----
// -------------------------------------------------------------------------

/// Last modified: 24.07.2013, P.B. 

#include <iostream>
#include "MpdConstField.h"
#include "MpdFieldMap.h"
#include "MpdFieldPar.h"
#include "FairParamList.h"
#include "MpdMultiField.h"
#include "TObjArray.h"

using namespace std;

// ------   Constructor   --------------------------------------------------
MpdFieldPar::MpdFieldPar(const char* name, const char* title,
			 const char* context) 
  : FairParGenericSet(name, title, context) {
  fType = -1;
  fXmin = fXmax = fYmin = fYmax = fZmin = fZmax = 0.;
  fBx   = fBy   = fBz   = 0.;
  fMapName = "";
  fPosX = fPosY = fPosZ = 0.;
  fScale = 0.;
}

// -------------------------------------------------------------------------
MpdFieldPar::MpdFieldPar() 
{
  fType = -1;
  fXmin = fXmax = fYmin = fYmax = fZmin = fZmax = 0.;
  fBx   = fBy   = fBz   = 0.;
  fMapName = "";
  fPosX = fPosY = fPosZ = 0.;
  fScale = 0.;
}

// ------   Destructor   ---------------------------------------------------
MpdFieldPar::~MpdFieldPar() { }

// ------   Put parameters   -----------------------------------------------
void MpdFieldPar::putParams(FairParamList* list) {

  if ( ! list ) return;
  
  list->add("Field Type", fType);

  if ( fType == 0 ) {                    // constant field
    list->add("Field min x", fXmin);
    list->add("Field max x", fXmax);
    list->add("Field min y", fYmin);
    list->add("Field max y", fYmax);
    list->add("Field min z", fZmin);
    list->add("Field max z", fZmax);
    list->add("Field Bx", fBx);
    list->add("Field By", fBy);
    list->add("Field Bz", fBz);
  }
  
  else if (fType ==1) {    // field map
    list->add("Field map name", fMapName);
    list->add("Field x position", fPosX);
    list->add("Field y position", fPosY);
    list->add("Field z position", fPosZ);
    list->add("Field scaling factor", fScale);
  }
}

// --------   Get parameters   ---------------------------------------------
Bool_t MpdFieldPar::getParams(FairParamList* list) {
  
  if ( ! list ) return kFALSE;
  
  if ( ! list->fill("Field Type", &fType) ) return kFALSE;
  
  if ( fType == 0 ) {                    // constant field
    if ( ! list->fill("Field min x", &fXmin) ) return kFALSE;
    if ( ! list->fill("Field max x", &fXmax) ) return kFALSE;
    if ( ! list->fill("Field min y", &fYmin) ) return kFALSE;
    if ( ! list->fill("Field max y", &fYmax) ) return kFALSE;
    if ( ! list->fill("Field min z", &fZmin) ) return kFALSE;
    if ( ! list->fill("Field max z", &fZmax) ) return kFALSE;
    if ( ! list->fill("Field Bx", &fBx) ) return kFALSE;
    if ( ! list->fill("Field By", &fBy) ) return kFALSE;
    if ( ! list->fill("Field Bz", &fBz) ) return kFALSE;
  }
  
  else if (fType ==1) {    // field map
    Text_t mapName[80];
    if ( ! list->fill("Field map name", mapName, 80) ) return kFALSE;
    fMapName = mapName;
    if ( ! list->fill("Field x position", &fPosX) )  return kFALSE;
    if ( ! list->fill("Field y position", &fPosY) )  return kFALSE;
    if ( ! list->fill("Field z position", &fPosZ) )  return kFALSE;
    if ( ! list->fill("Field scaling factor", &fScale) ) return kFALSE;
  }
  return kTRUE;
}

// ---------   Set parameters from FairField   ------------------------------
void MpdFieldPar::SetParameters(FairField* field) {
  
  if ( ! field ) {
    cerr << "-W- MpdFieldPar::SetParameters: Empty field pointer!" << endl;
    return;
  }
  
  fType = field->GetType();
  
  if ( fType == 0 ) {                                 // constant field
    MpdConstField* fieldConst = (MpdConstField*) field;
    fBx = fieldConst->GetBx();
    fBy = fieldConst->GetBy();
    fBz = fieldConst->GetBz();
    fXmin = fieldConst->GetXmin();
    fXmax = fieldConst->GetXmax();
    fYmin = fieldConst->GetYmin();
    fYmax = fieldConst->GetYmax();
    fZmin = fieldConst->GetZmin();
    fZmax = fieldConst->GetZmax();
    fMapName = "";
    fPosX = fPosY = fPosZ = fScale = 0.;
  }
  
  else if ( fType ==1 ) {              // field map
    MpdFieldMap* fieldMap = (MpdFieldMap*) field;
    fBx = fBy = fBz = 0.;
    fXmin = fXmax = fYmin = fYmax = fZmin = fZmax = 0.;
    
    fMapName = field->GetName();
    fPosX   = fieldMap->GetPositionX();
    fPosY   = fieldMap->GetPositionY();
    fPosZ   = fieldMap->GetPositionZ();
    fScale  = fieldMap->GetScale();
  } 
  
  else if (fType == 5) {
    MpdMultiField *fMulti=(MpdMultiField *)field;
    TObjArray *fieldlist = fMulti->GetFieldList();
    TIterator* FieldIter = fieldlist->MakeIterator();
    FieldIter->Reset();
    FairField *fi=0;
    while( (fi = (FairField*)FieldIter->Next() ) ) {
      SetParameters(fi);
    }
    
    delete  FieldIter;
  }
  else {
    cerr << "-W- MpdFieldPar::SetParameters: Unknown field type "
	 << fType << "!" << endl;
    fBx = fBy = fBz = 0.;
    fXmin = fXmax = fYmin = fYmax = fZmin = fZmax = 0.;
    fMapName = "";
    fPosX = fPosY = fPosZ = fScale = 0.;
  }
  
  return;

}

ClassImp(MpdFieldPar)

