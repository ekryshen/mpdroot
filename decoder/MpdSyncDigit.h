//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPDSyncDIGIT_H
#define __MPDSyncDIGIT_H

#include <iostream>

#include "TObject.h"
//------------------------------------------------------------------------------------------------------------------------
class MpdSyncDigit : public TObject 
{
    UInt_t fSerial = 0;
    Long64_t fGlobalEvent = 0;
    Long64_t fT_sec = 0;
    Long64_t fT_ns = 0;

public:
    MpdSyncDigit() = default;
    MpdSyncDigit(UInt_t iSerial, Long64_t iEvent, Long64_t t_sec, Long64_t t_ns);
    virtual ~MpdSyncDigit() = default;

    void 		Dump() const;
    UInt_t GetSerial() const { return fSerial; }
    Long64_t GetEvent() const { return fGlobalEvent; }
    Long64_t GetTime_sec()const { return fT_sec; }
    Long64_t GetTime_ns() const { return fT_ns; }

    ClassDef(MpdSyncDigit, 1);
};
#endif 
//------------------------------------------------------------------------------------------------------------------------

