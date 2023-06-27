//------------------------------------------------------------------------------------------------------------------------
/// \class MpdRawDataDecoder
/// 
/// \brief 
/// \author Victor Baryshnikov
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------
#include "MpdEventHeader.h"

ClassImp(MpdEventHeader)

//------------------------------------------------------------------------------------------------------------------------
MpdEventHeader::MpdEventHeader(UInt_t period, UInt_t run,  UInt_t rack, UInt_t event)
: m_PeriodId(period), m_RackId(rack), m_EventId(event)
{   
 	SetRunId(run);
}
//------------------------------------------------------------------------------------------------------------------------
void MpdEventHeader::Print(const char *comment, std::ostream &os) 
{
   if (comment != nullptr) os << comment;

   os<<" [MpdEventHeader] PeriodId="<<printNan(m_PeriodId)<<", runId="<<GetRunId()<<", rackId="<<printNan(m_RackId)<<", eventId="<<printNan(m_EventId)
	<<", time="<<GetEventTime()<<", fileId="<<GetInputFileId();
}
//------------------------------------------------------------------------------------------------------------------------

