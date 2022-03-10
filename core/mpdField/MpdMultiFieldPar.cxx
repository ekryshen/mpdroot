// -------------------------------------------------------------------------
//                            MpdMultiFieldPar source file             -----
//        Created 23/07/13  by P. Batyuk (MPD) <batyuk@jinr.ru>        -----
//                          from MpdMultiFieldPar (PNDROOT)            -----
// -------------------------------------------------------------------------

/// Last modified: 24.07.2013, P.B.

#include "MpdMultiFieldPar.h"
#include "MpdMultiFieldPar.h"
#include "MpdConstField.h"
#include "MpdConstPar.h"
#include "FairParamList.h"
#include "MpdMapPar.h"
#include "FairRuntimeDb.h"
#include "FairRun.h"
#include "MpdFieldMap.h"
#include "MpdFieldPar.h"

// ------   Constructor   --------------------------------------------------
MpdMultiFieldPar::MpdMultiFieldPar(const char* name, const char* title, const char* context) 
  : MpdMapPar(name, title, context) 
{
  
  fParArray=new TObjArray();
}

//-----------------------------------------------------------------------
MpdMultiFieldPar::MpdMultiFieldPar() { }

//-----------------------------------------------------------------------
MpdMultiFieldPar::~MpdMultiFieldPar() { }

//-----------------------------------------------------------------------
void MpdMultiFieldPar::putParams(FairParamList* list)
{
  if ( ! list ) return;
  list->addObject("List of Field par", fParArray);
  list->add("Field Type", fType);
}

//-----------------------------------------------------------------------
Bool_t MpdMultiFieldPar::getParams(FairParamList* l)
{
  if (!l->fillObject("list of fields Par",fParArray))return kFALSE;
  if ( ! l->fill("Field Type", &fType) ) return kFALSE;
  
  return kTRUE;
} 

//-----------------------------------------------------------------------
void MpdMultiFieldPar:: SetParameters(FairField* field)
{	 
  fType=5;
  FairRuntimeDb *rtdb=FairRuntimeDb::instance();
  FairRun *fRun= FairRun::Instance();
  MpdMultiField *fmField = (MpdMultiField *)field;
  TObjArray *fArray=fmField->GetFieldList();
  TIterator *Iter=fArray->MakeIterator();
  Iter->Reset();
  FairField* fField = NULL;
  Int_t Type=-1;
  while( (fField = (FairField*)Iter->Next() ) ) {
    Type=fField->GetType();
    if(Type==0) {
      MpdConstField *fc= (MpdConstField *)fField;
      MpdConstPar *cp = (MpdConstPar*) rtdb->getContainer("MpdConstPar");
      cp->SetParameters(fc);
      cp->setInputVersion(fRun->GetRunId(),1);
      fParArray->AddLast(cp);
    }
    
    else if(Type==1) {
      MpdFieldMap *fm= (MpdFieldMap*)fField;
      MpdFieldPar* fmp = (MpdFieldPar*) rtdb->getContainer("MpdFieldPar");
      fmp->SetParameters(fm);
      fmp->setInputVersion(fRun->GetRunId(),1);
      fParArray->AddLast(fmp);
    } 
  }
  delete  Iter;
}                      

ClassImp(MpdMultiFieldPar)
  
