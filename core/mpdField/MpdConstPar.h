// -------------------------------------------------------------------------
//                            MpdContPar header file                   -----
//        Created 23/07/13  by P. Batyuk (MPD) <batyuk@jinr.ru>        -----
//                          from MpdContPar (PNDROOT)                  -----
// -------------------------------------------------------------------------

#ifndef MPDCONSTPAR_H
#define MPDCONSTPAR_H 1
#include "MpdMapPar.h"

class FairParamList;

class MpdConstPar : public MpdMapPar
{
  
 public:
  
  /** Standard constructor  **/
  MpdConstPar(const char* name, const char* title, const char* context);
  
  /** default constructor  **/
  MpdConstPar();
  
  /** Destructor **/
  ~MpdConstPar();
  
  void putParams(FairParamList* list);

  /** Get parameters **/
  Bool_t getParams(FairParamList* list);
  
  /** Set parameters from FairField  **/
  void SetParameters(FairField* field);
  
  Double_t GetBx()        const { return fBx; }
  Double_t GetBy()        const { return fBy; }
  Double_t GetBz()        const { return fBz; }

 protected:
  
  /** Field values in [kG] **/
  Double_t fBx, fBy, fBz;
  
  ClassDef(MpdConstPar,1);
};

#endif
