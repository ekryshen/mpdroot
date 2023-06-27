//------------------------------------------------------------------------------------------------------------------------
/// \class MpdEventHeader
/// 
/// \brief 
/// \author Victor Baryshnikov
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_EVENT_HEADER_H
#define __MPD_EVENT_HEADER_H 

#include <iostream>

#include "FairRootManager.h"
#include "FairEventHeader.h"
//------------------------------------------------------------------------------------------------------------------------
class MpdEventHeader : public FairEventHeader 
{

class printNan
{
UInt_t value;
public:
	printNan(UInt_t v): value(v){};
	friend std::ostream& operator<<(std::ostream& os, const printNan& dt){ if(dt.value == -1) os<<"NaN"; else os<<dt.value; return os;}
};

public:
	UInt_t 	m_PeriodId = (UInt_t) -1;
	UInt_t 	m_RackId = (UInt_t) -1;
	UInt_t 	m_EventId = (UInt_t) -1;

  	MpdEventHeader() = default;
    	MpdEventHeader(UInt_t period, UInt_t run, UInt_t rack, UInt_t event);
    	virtual ~MpdEventHeader() = default;
   
   	virtual void Register(Bool_t Persistence = kTRUE)
   	{
        	FairRootManager::Instance()->Register("MpdEventHeader.", "EvtHeader", this, Persistence);
   	}    

	// copy assignment
	MpdEventHeader& operator=(const MpdEventHeader& other) = default;

    	void 	Print(const char *comment = nullptr, std::ostream &os = std::cout);

ClassDef(MpdEventHeader, 1)
};

#endif 
//------------------------------------------------------------------------------------------------------------------------

