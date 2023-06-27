//------------------------------------------------------------------------------------------------------------------------
#include "TLVBlock.h"
#include "MpdTTVXSDigit.h"

ClassImp(MpdTTVXSDigit)
//------------------------------------------------------------------------------------------------------------------------
MpdTTVXSDigit::MpdTTVXSDigit()
: fSerial(0), fTriggerType(0), fTriggerSource(0)
{ 
    fEventTime.SetSec(0);
    fEventTime.SetNanoSec(0);   
}
//------------------------------------------------------------------------------------------------------------------------
MpdTTVXSDigit::MpdTTVXSDigit(UInt_t serial, TTimeStamp time, UShort_t type, UShort_t source)
: fSerial(serial), fEventTime(time), fTriggerType(type), fTriggerSource(source)
{

}
//------------------------------------------------------------------------------------------------------------------------
void 		MpdTTVXSDigit::Dump() const
{
	std::cout<<"\n MpdTTVXSDigit: fSerial="<<toHex(fSerial, 4)<<", fEventTime="<<fEventTime.AsDouble()<<", fTriggerType="<<fTriggerType<<", fTriggerSource="<<fTriggerSource;
}
//------------------------------------------------------------------------------------------------------------------------    
        
