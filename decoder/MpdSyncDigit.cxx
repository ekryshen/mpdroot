//------------------------------------------------------------------------------------------------------------------------
#include "TLVBlock.h"
#include "MpdSyncDigit.h"

ClassImp(MpdSyncDigit)
//------------------------------------------------------------------------------------------------------------------------
MpdSyncDigit::MpdSyncDigit(UInt_t serial, Long64_t event, Long64_t sec, Long64_t ns) 
: fSerial(serial), fGlobalEvent(event), fT_sec(sec), fT_ns(ns)
{
    
}
//------------------------------------------------------------------------------------------------------------------------
void 		MpdSyncDigit::Dump() const
{
	std::cout<< " MpdSyncDigit: fSerial="<<toHex(fSerial, 4)<<", fGlobalEvent="<<fGlobalEvent<<", fT_sec= "<<fT_sec <<", fT_ns= "<<fT_ns<<std::endl;
}
//------------------------------------------------------------------------------------------------------------------------

