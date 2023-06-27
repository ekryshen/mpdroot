//------------------------------------------------------------------------------------------------------------------------
/// \class MpdTofDigitProducer
/// 
/// \brief 
/// \author Victor Baryshnikov
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------

#include <TClonesArray.h>

#include "FairRunAna.h"
#include "MpdEventHeader.h"
#include "MpdTofDigitImpl.h"

#include "MpdTofDigitProducer.h"

using namespace std;

ClassImp(MpdTofDigitProducer)
//------------------------------------------------------------------------------------------------------------------------
MpdTofDigitProducer::MpdTofDigitProducer(FairRunAna* run, TOF_DIGIT_PRODUCER_DESC& desc, Int_t verbose, const char *name)
  : FairTask(name, verbose), pImpl(std::make_unique<Impl>(this, verbose)), m_Run(run)
{
	m_desc = desc;
	
}
//------------------------------------------------------------------------------------------------------------------------
MpdTofDigitProducer::~MpdTofDigitProducer()
{

}
//------------------------------------------------------------------------------------------------------------------------
InitStatus MpdTofDigitProducer::Init()
{
	if(m_flnm.empty()) // input device digits from FairRunAna chain.
	{
      		aSync 	= (TClonesArray*) FairRootManager::Instance()->GetObject("SyncDigit");
      		aTdc 	= (TClonesArray*) FairRootManager::Instance()->GetObject("TDCDigit");
      		aTtvxs 	= (TClonesArray*) FairRootManager::Instance()->GetObject("TTVXSDigit");

	}
assert(aSync);
assert(aTdc);
assert(aTtvxs);

        // Create and register output array
        aTofDigits = new TClonesArray("MpdTofDigit");
	FairRootManager::Instance()->Register("TOFDigit", "Tof", aTofDigits, kTRUE);

	auto status = kSUCCESS;
	
	// Initialize Tof digitizer implementation
	if(m_desc.rackId != -1 && m_desc.runId != -1)
	{
 		status = pImpl->Init(m_desc);
		if(status == kSUCCESS) LOG(info)<<"[MpdTofDigitProducer::Init] Initialization finished succesfully.";
	}
	
return status;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofDigitProducer::Exec(Option_t *option)
{
	// Cleanup containers
	aTofDigits->Clear();
	
	auto header = dynamic_cast<MpdEventHeader*>(m_Run->GetEventHeader());
	if(header)
	{
		// copy periodId from config to MpdEventHeader
		header->m_PeriodId = m_desc.periodId;

		// delay initialization by first time
		if(!pImpl->IsReady())
		{	
			// save real runId and rackId from MpdEventHeader to config
			m_desc.runId = header->GetRunId();
			m_desc.rackId = header->m_RackId;

			// MpdTofDigitProducer::Impl::Init MUST be AFTER MpdRawDataDecoder::Exec(read real RunId and rackId from file) ---> delay init
			auto status = pImpl->Init(m_desc);

assert(status == kSUCCESS);
		}

		pImpl->Exec(m_desc);
	
		LOG(debug1)<<"[MpdTofDigitProducer::Exec] tof digits="<<aTofDigits->GetEntries()<<", aTdc digits="<<aTdc->GetEntries();
	}
	else LOG(error)<<"MpdEventHeader don't exist";
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofDigitProducer::Finish()
{

}
//------------------------------------------------------------------------------------------------------------------------

