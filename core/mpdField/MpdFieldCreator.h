// -------------------------------------------------------------------------
//                            MpdFieldCreator header file              -----
//        Created 23/07/13  by P. Batyuk (MPD) <batyuk@jinr.ru>        -----
//                          from MpdFieldCreator (PNDROOT)             -----
// -------------------------------------------------------------------------

#ifndef MPDFIELDCREATOR_H
#define MPDFIELDCREATOR_H

#include "FairField.h"
#include "FairFieldFactory.h"
#include "MpdFieldPar.h"

class  MpdSolenoidPar;
class  MpdDipolePar  ;
class  MpdTransPar   ;
class  MpdConstPar   ;
class  MpdMultiFieldPar ;

class MpdFieldCreator : public FairFieldFactory
{
 
 public:
  MpdFieldCreator();
  virtual ~MpdFieldCreator();
  virtual FairField* createFairField();
  virtual void SetParm();
  ClassDef(MpdFieldCreator,1);

 protected:
   
  MpdFieldPar* fFieldPar;
  MpdSolenoidPar    *fSPar;
  MpdDipolePar      *fDPar;
  MpdTransPar       *fTPar;
  MpdConstPar       *fCPar;
  MpdMultiFieldPar  *fMPar;
};
#endif 
