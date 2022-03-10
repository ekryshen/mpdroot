// -------------------------------------------------------------------------
//                            MpdFieldCreator source file              -----
//        Created 23/07/13  by P. Batyuk (MPD) <batyuk@jinr.ru>        -----
//                          from MpdFieldCreator (PNDROOT)             -----
// -------------------------------------------------------------------------

/// Last modified: 24.07.2013, P.B.

#include "MpdFieldCreator.h"
#include "MpdFieldMap.h"
#include "MpdConstField.h"
#include "MpdMultiField.h"
#include "MpdMultiFieldPar.h"
#include "MpdConstPar.h"
#include "MpdMapPar.h"

#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "TObjArray.h"
#include <iostream>
using std::cout;
using std::endl;

static MpdFieldCreator gMpdFieldCreator;

//-------------------------------------------------------------------
MpdFieldCreator::MpdFieldCreator()
        :FairFieldFactory()
{
	fCreator=this;
  fSPar=0;
  fDPar=0;
  fTPar=0;
  fCPar=0;
  fMPar=0;
  fFieldPar=0;
}

//-------------------------------------------------------------------
MpdFieldCreator::~MpdFieldCreator() { }

//-------------------------------------------------------------------
void MpdFieldCreator::SetParm()
{
  cout << "MpdFieldCreator::SetParm()" <<endl;
  FairRunAna *Run = FairRunAna::Instance();
  FairRuntimeDb *RunDB = Run->GetRuntimeDb();
  
  fFieldPar = (MpdFieldPar*) RunDB->getContainer("MpdFieldPar");
  fCPar = (MpdConstPar     *)RunDB->getContainer("MpdConstPar");
  fMPar = (MpdMultiFieldPar*)RunDB->getContainer("MpdMultiFieldPar");
}

//-------------------------------------------------------------------
FairField* MpdFieldCreator::createFairField()
{       
  cout << "MpdFieldCreator::createFairField()" <<endl;
  FairField *fMagneticField=0;
  MpdMultiField *MField=0;
  Int_t Type=-1;
  Bool_t multi=kFALSE;
  MField = new MpdMultiField(fMPar);
  Type= fMPar->GetType();
  if (Type==-1) {delete MField; MField=0;}
  if(MField) {
    multi=kTRUE;
    TObjArray *fParArray=fMPar->GetParArray();
    TIterator *Iter=fParArray->MakeIterator();
    Iter->Reset();
    MpdMapPar* fPar = NULL;
    while( (fPar = (MpdMapPar*)Iter->Next() ) ) {
      fPar->Print();
 
      if (fPar->GetType()==0) {
	FairField *fField1 = new MpdConstField((MpdConstPar*)fPar);
	MField->AddField(fField1);
      }
      
      else if (fPar->GetType()==1) {
	FairField *fField2 = new MpdFieldMap((MpdFieldPar*)fPar);
	MField->AddField(fField2);
      }
    }
    delete Iter;
  }
  
  FairField *fField1 = new MpdConstField(fCPar);
  Type= fCPar->GetType();
  if (Type==-1) {delete fField1; fField1=0;}
  if(fField1)fMagneticField=fField1;	
  
  FairField *fField2 = new MpdFieldMap(fFieldPar);
  Type= fFieldPar->GetType();
  if (Type==-1){delete fField2; fField2=0;}
  if(fField2)fMagneticField=fField2;	
  
  if (multi) {
    MField->Init();
    return MField;
  } 
  else{
    if(fMagneticField){
      fMagneticField->Init();	 
    }
    return fMagneticField;
  }
}

ClassImp(MpdFieldCreator)

