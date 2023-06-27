//------------------------------------------------------------------------------------------------------------------------
#include "TLVBlock.h"
#include "MpdTDCDigit.h"

ClassImp(MpdTDCDigit)
//------------------------------------------------------------------------------------------------------------------------
MpdTDCDigit::MpdTDCDigit(UInt_t serial, UChar_t type, Bool_t leading, UChar_t channel, UInt_t value) 
: fSerial(serial), fType(type), fLeading(leading), fChannel(channel), fValue(value)
{

}
//------------------------------------------------------------------------------------------------------------------------
void 		MpdTDCDigit::Dump() const
{
	std::cout<<" MpdTDCDigit: fSerial="<<toHex(fSerial, 4)<<", fType="<<toHex((int)fType, 2) 
		<< ", fValue=" << fValue<<", fChannel= "<<(int)fChannel <<", fLeading= "<<fLeading<<std::endl;
}
//------------------------------------------------------------------------------------------------------------------------
