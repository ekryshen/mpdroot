// -------------------------------------------------------------------------
//                            MpdMultiField source file                -----
//        Created 23/07/13  by P. Batyuk (MPD) <batyuk@jinr.ru>        -----
//                          from MpdMultiField    (PNDROOT)            -----
// -------------------------------------------------------------------------

/// Last modified: 24.07.2013, P.B.

#include <iomanip>
#include <iostream>
#include <fstream>

#include "MpdMultiField.h"
#include "MpdRegion.h"
#include "MpdConstField.h"
#include "MpdFieldMap.h"
#include "TObjArray.h"
#include "MpdMapPar.h"
#include "MpdMultiFieldPar.h"


using namespace std;

// -------------   Default constructor  ----------------------------------
MpdMultiField::MpdMultiField() {

    fMaps= new TObjArray(10);
    fNoOfMaps=0;
    fType = 5;
}

// ------------   Constructor from MpdFieldPar   --------------------------
MpdMultiField::MpdMultiField(MpdMultiFieldPar* fieldPar) {
   fType = 5;
   fMaps= new TObjArray(10);
   fNoOfMaps=0;
   TObjArray *fArray= fieldPar->GetParArray();
   if(fArray->IsEmpty()) fType=-1;
}

// ------------   Destructor   --------------------------------------------
MpdMultiField::~MpdMultiField() { }

// -----------   Adding fields   ------------------------------------------
void MpdMultiField::AddField(FairField *field) {
	
   if(field){
      fMaps->AddLast(field);	
      fNoOfMaps++;
   }
}

// -----------   Intialisation   ------------------------------------------
void MpdMultiField::Init() {
  MpdConstField *field=0;
  MpdFieldMap *fieldMap=0;
  for (Int_t n=0; n<=fNoOfMaps; n++){
    fieldMap = dynamic_cast<MpdFieldMap *>(fMaps->At(n));
    field = dynamic_cast<MpdConstField *>(fMaps->At(n));
    if(fieldMap){
      fieldMap->Init();
      fFieldMaps.insert( pair<MpdRegion*, FairField*>(new MpdRegion(fieldMap->GetZmin(),fieldMap->GetZmax()) , fieldMap ));
    }else if(field){
      field->Init();
      fFieldMaps.insert( pair<MpdRegion*, FairField*>(new MpdRegion(field->GetZmin(),field->GetZmax() ), field ));
    }
  }
}

// -----------   Get x component of the field   ---------------------------
Double_t MpdMultiField::GetBx(Double_t x, Double_t y, Double_t z) {
  MpdRegion *fReg=0;
  FairField *fField=0;
  std::map <MpdRegion*, FairField* >::iterator fMapIter;
  for (fMapIter=fFieldMaps.begin(); fMapIter!= fFieldMaps.end();fMapIter++ ){
    fReg=fMapIter->first;
    if(fReg->IsInside(z)){
      fField=fMapIter->second;
      break;
    }
  }
  if(fField)  return  fField->GetBx( x, y,  z);	
  else return 0;
}

// -----------   Get y component of the field   ---------------------------
Double_t MpdMultiField::GetBy(Double_t x, Double_t y, Double_t z) {
  
  MpdRegion *fReg=0;
  FairField *fField=0;
  std::map <MpdRegion*, FairField* >::iterator fMapIter;
  for (fMapIter=fFieldMaps.begin(); fMapIter!= fFieldMaps.end();fMapIter++ ){
    fReg=fMapIter->first;
    if(fReg->IsInside(z)){
      fField=fMapIter->second;
      break;
    }
  }
  if(fField)  return fField->GetBy( x, y,  z);
  else return 0;
}

// -----------   Get z component of the field   ---------------------------
Double_t MpdMultiField::GetBz(Double_t x, Double_t y, Double_t z) {
  
  MpdRegion *fReg=0;
  FairField *fField=0;
  std::map <MpdRegion*, FairField* >::iterator fMapIter;
  for (fMapIter=fFieldMaps.begin(); fMapIter!= fFieldMaps.end();fMapIter++ ){
    fReg=fMapIter->first;
    if(fReg->IsInside(z)){
      fField=fMapIter->second;
      break;
    }
  }
  if(fField) return fField->GetBz( x, y,  z);
  else return 0;
}

//---------------------------------------------------------------------------
void MpdMultiField::GetFieldValue(const Double_t point[3], Double_t* bField)
{
  
  MpdRegion *fReg=0;
  FairField *fField=0;
  std::map <MpdRegion*, FairField* >::iterator fMapIter;
  for (fMapIter=fFieldMaps.begin(); fMapIter!= fFieldMaps.end();fMapIter++ ){
    fReg=fMapIter->first;
    if(fReg->IsInside(point[2])){
      fField=fMapIter->second;
      break;
    }
  }
  if(fField){
    /* bField[0] = fField->GetBx(point[0], point[1], point[2]);
       bField[1] = fField->GetBy(point[0], point[1], point[2]);
       bField[2] = fField->GetBz(point[0], point[1], point[2]);
    */
    fField->GetBxyz(point, bField);
  }
  else{
    bField[0] = 0;
    bField[1] = 0;
    bField[2] = 0;
  }
}

// ---------   Screen output   --------------------------------------------
void MpdMultiField::Print() {  
  for (Int_t n=0; n<=fNoOfMaps; n++){
    FairField *fieldMap = dynamic_cast<FairField *>(fMaps->At(n));
    if(fieldMap) fieldMap->Print("");
  }
}

ClassImp(MpdMultiField)
  
