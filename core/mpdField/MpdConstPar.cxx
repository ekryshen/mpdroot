// -------------------------------------------------------------------------
//                            MpdContPar source file                   -----
//        Created 23/07/13  by P. Batyuk (MPD) <batyuk@jinr.ru>        -----
//                          from MpdContPar (PNDROOT)                  -----
// -------------------------------------------------------------------------

/// Last modified: 24.07.2013, P.B.

#include "MpdConstPar.h"
#include "MpdConstField.h"
#include "FairParamList.h"

// ------   Constructor   --------------------------------------------------
MpdConstPar::MpdConstPar(const char* name, const char* title, const char* context) 
  : MpdMapPar(name, title, context) 
{
   fBx   = fBy   = fBz   = 0.;
   fType=0;
}

//-----------------------------------------------------------------
MpdConstPar::MpdConstPar() 
{
   fBx   = fBy   = fBz   = 0.;
   fType=0;
}

//-----------------------------------------------------------------
MpdConstPar::~MpdConstPar() { }

// ----------------------------------------------------------------
void MpdConstPar::putParams(FairParamList* list)
{
  if ( ! list ) return;

  list->add("Field Type", fType);
 
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

//-----------------------------------------------------------------
Bool_t MpdConstPar::getParams(FairParamList* list)
{
  if ( ! list ) return kFALSE;
  if ( ! list->fill("Field Type", &fType) ) return kFALSE;
  if ( ! list->fill("Field min x", &fXmin) ) return kFALSE;
  if ( ! list->fill("Field max x", &fXmax) ) return kFALSE;
  if ( ! list->fill("Field min y", &fYmin) ) return kFALSE;
  if ( ! list->fill("Field max y", &fYmax) ) return kFALSE;
  if ( ! list->fill("Field min z", &fZmin) ) return kFALSE;
  if ( ! list->fill("Field max z", &fZmax) ) return kFALSE;
  if ( ! list->fill("Field Bx", &fBx) ) return kFALSE;
  if ( ! list->fill("Field By", &fBy) ) return kFALSE;
  if ( ! list->fill("Field Bz", &fBz) ) return kFALSE;
  return kTRUE;
} 

//------------------------------------------------------------------
void MpdConstPar:: SetParameters(FairField* field)
{
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
    fType = fieldConst->GetType();
}

ClassImp(MpdConstPar)
