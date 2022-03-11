//
// Holds information about hit in FHCal
//

#include "MpdMiniMessMgr.h"
#include "MpdMiniFHCalHit.h"

ClassImp(MpdMiniFHCalHit)

   //_________________
   MpdMiniFHCalHit::MpdMiniFHCalHit()
   : TObject(), fId(0), fEDep(0) /*,fAdc(0)*/
{
   // Default constructor
}

//_________________
MpdMiniFHCalHit::MpdMiniFHCalHit(const MpdMiniFHCalHit &hit)
   : TObject(), fId(hit.fId), fEDep(hit.fEDep) /*, fAdc(hit.fAdc) */
{
   // Copy constructor
}

//_________________
MpdMiniFHCalHit::~MpdMiniFHCalHit()
{
   // Empty destructor
}

//_________________
void MpdMiniFHCalHit::Print(const Char_t *option __attribute__((unused))) const
{
   LOG_INFO << "DetId: " << detector() << " ModId: " << module() << " ChanId: " << channel()
            << " eDep: " << eDep() /*<<" adc: " << adc() */
            << endm;
}
