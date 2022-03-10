// -------------------------------------------------------------------------
//                            MpdFieldContFact header file             -----
//        Created 23/07/13  by P. Batyuk (MPD) <batyuk@jinr.ru>        -----
//                          from MpdFieldContFact (PNDROOT)            -----
// -------------------------------------------------------------------------

#ifndef MPDFIELDCONTFACT_H
#define MPDFIELDCONTFACT_H 1

#include "FairContFact.h"

class FairContainer;
class FairParSet;

class MpdFieldContFact : public FairContFact
{
  
 public:
  
  /** Constructor **/
  MpdFieldContFact();

  /** Destructor **/
  ~MpdFieldContFact();
  
  /** Create containers
   ** Creates the requested parameter sets (MpdFieldPar) 
   **/
  FairParSet* createContainer(FairContainer* container);
  
 private:
  
  /** Set all containers  
   ** Creates container objects with all accepted contexts and adds them
   ** to the list of containers for the field library. 
   **/
  void SetAllContainers();
  
  ClassDef(MpdFieldContFact,1);
};

#endif
