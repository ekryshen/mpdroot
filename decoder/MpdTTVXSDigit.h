//------------------------------------------------------------------------------------------------------------------------
/* 
 * File:   MpdTTVXSDigit.h
 * Author: mikhailr
 *
 * Created on September 15, 2021, 11:09 AM
 */

#ifndef __MPDTTVXSDIGIT_H
#define __MPDTTVXSDIGIT_H

#include <iostream>

#include "TTimeStamp.h"
#include "TObject.h"
//------------------------------------------------------------------------------------------------------------------------
class MpdTTVXSDigit : public TObject 
{
    UInt_t 	fSerial;
    TTimeStamp 	fEventTime;
    UShort_t 	fTriggerType;
    UShort_t 	fTriggerSource;

public:
    MpdTTVXSDigit();
    MpdTTVXSDigit(UInt_t serial, TTimeStamp time, UShort_t type, UShort_t source);
    virtual ~MpdTTVXSDigit() = default;

    void 		Dump() const;
    UInt_t 	GetSerial() const { return fSerial;  }
    TTimeStamp 	GetTime() const { return fEventTime; }
    UShort_t 	GetTriggerType() const { return fTriggerType; }
    UShort_t 	GetTriggerSourse() const { return fTriggerSource; }

    ClassDef(MpdTTVXSDigit, 1)
};
#endif 
//------------------------------------------------------------------------------------------------------------------------

