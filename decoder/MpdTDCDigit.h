//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPDTDCDIGIT_H
#define __MPDTDCDIGIT_H 1

#include <iostream>

#include "TObject.h"
//------------------------------------------------------------------------------------------------------------------------
class MpdTDCDigit : public TObject 
{
    UInt_t 	fSerial = 0;
    UInt_t 	fValue = 0;
    UChar_t 	fType = 0;
    UChar_t 	fChannel  = 0;
    Bool_t 	fLeading = 0;

public:
    MpdTDCDigit() = default;
    MpdTDCDigit(UInt_t iSerial, UChar_t iType, Bool_t iLeading, UChar_t iChannel, UInt_t iValue);
    virtual ~MpdTDCDigit() = default;

    void 		Dump() const;
    UInt_t 	GetSerial() const { return fSerial; }
    UChar_t 	GetType() const { return fType; }
    Bool_t 	GetLeading() const { return fLeading; }
    UChar_t 	GetChannel() const { return fChannel; }
    UInt_t 	GetValue() const { return fValue; }

    bool operator< ( const MpdTDCDigit& rhs)const // sort by guid(serial + channel) 
    {
	if(this->fSerial != rhs.fSerial) return (this->fSerial < rhs.fSerial);	
	
	return (this->fChannel < rhs.fChannel);					
    }

    ClassDef(MpdTDCDigit, 1);
};
#endif 
//------------------------------------------------------------------------------------------------------------------------

