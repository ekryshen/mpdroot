// -------------------------------------------------------------------------
//                            MpdMultiFieldPar header file             -----
//        Created 23/07/13  by P. Batyuk (MPD) <batyuk@jinr.ru>        -----
//                          from MpdMultiFieldPar (PNDROOT)            -----
// -------------------------------------------------------------------------

#ifndef MPDMULTIFIELDPAR_H
#define MPDMULTIFIELDPAR_H 1

#include "MpdMapPar.h"
#include "MpdMultiField.h"

class FairParamList;

class MpdMultiFieldPar : public MpdMapPar
{
  
 public:
  
  /** Standard constructor  **/
  MpdMultiFieldPar(const char* name, const char* title, const char* context);
  
  /** default constructor  **/
  MpdMultiFieldPar();
  
  /** Destructor **/
  ~MpdMultiFieldPar();
  
  void putParams(FairParamList* list);
  
  /** Get parameters **/
  Bool_t getParams(FairParamList* list);
  
  /** Set parameters from FairField  **/
  void SetParameters(FairField* field);
  
  TObjArray *GetParArray(){return fParArray; }
  
 protected:
  
  TObjArray *fParArray;
  
  ClassDef(MpdMultiFieldPar,1);
  
};

#endif
