//------------------------------------------------------------------------------------------------------------------------
/// \class MpdRawDataDecoder
/// 
/// \brief 
/// \author Victor Baryshnikov
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------
#include <assert.h>

#include <TClonesArray.h>

#include "FairLogger.h"
#include "FairRunAna.h"
#include "FairEventHeader.h"
#include "MpdEventHeader.h"

// used digits:
#include "MpdSyncDigit.h"
#include "MpdTDCDigit.h"
#include "MpdTTVXSDigit.h"

// used devices:
#include "TDC72VXS.h"
#include "TTVXS.h"

#include "MpdRawDataDecoder.h"
#include "MpdDecoderImpl.h"

using namespace std;

ClassImp(MpdRawDataDecoder)
//------------------------------------------------------------------------------------------------------------------------
MpdRawDataDecoder::MpdRawDataDecoder(FairRunAna* run, const char *name, Int_t verbose, SAVE_BRANCH_FLAGS flags, RAW_DATA_FORMATS type)
  : FairTask(name, verbose), pImpl(std::make_unique<Impl>(this, type, verbose)), fSaveFlags(flags), m_Run(run)
{
	
}
//------------------------------------------------------------------------------------------------------------------------
MpdRawDataDecoder::~MpdRawDataDecoder()
{
	
}
//------------------------------------------------------------------------------------------------------------------------
InitStatus 	MpdRawDataDecoder::Init()
{
	// Create and register output arrays
	aSync = make_shared<TClonesArray>("MpdSyncDigit");
	aTdc = make_shared<TClonesArray>("MpdTDCDigit");
	aTtvxs = make_shared<TClonesArray>("MpdTTVXSDigit");

	FairRootManager::Instance()->Register("SyncDigit", "MpdDigits", aSync.get(), fSaveFlags & SAVE_SYNC_DIGITS);
	FairRootManager::Instance()->Register("TDCDigit", "MpdDigits", aTdc.get(), fSaveFlags & SAVE_TDC_DIGITS);
	FairRootManager::Instance()->Register("TTVXSDigit", "MpdDigits", aTtvxs.get(), fSaveFlags & SAVE_TTVXS_DIGITS);

	// Initialize raw data decoder implementation
	auto device = pImpl->CreateDevice(std::make_shared<TDC72VXS>(aTdc, aSync, fVerbose)); device->SetTitle("title");
	pImpl->CreateDevice(std::make_shared<TTVXS>(aTtvxs, fVerbose));

	auto status = pImpl->Init();
	LOG(info)<<"[MpdRawDataDecoder::Init] Initialization finished succesfully.";

return status;
}
//------------------------------------------------------------------------------------------------------------------------
void 			MpdRawDataDecoder::Exec(Option_t *option)
{
	// cleanup containers
   	aSync->Clear();
   	aTdc->Clear();
   	aTtvxs->Clear();

	// fill containers by new event data
	size_t eventId, position;
	auto status = pImpl->ReadNextEvent(eventId, position);

	if(status.first == kOk)
	{
/*		if(fVerbose > 3) 
		{
			DumpArray(aSync.get());
			DumpArray(aTdc.get());
			DumpArray(aTtvxs.get());
			cout<<endl;
		}
*/
	}
	else
	{
		LOG(debug)<<"pImpl->ReadNextEvent() return error status=("<<status.first<<", "<<status.second<<"), eventId="<<eventId<<", offset="<<position;
	}

	// Stop run
	if(status.second == kFinished)
	{
	 	m_Run->GetMainTask()->FinishTask();
		FinishTask();
		m_Run->StopProcessingLMD();
	}

	LOG(debug1)<<"[MpdRawDataDecoder::Exec] eventId="<<eventId<<", aSync= "<<aSync->GetEntries()<<", aTdc= "<<aTdc->GetEntries()<<", aTtvxs= "<<aTtvxs->GetEntries();
}
//------------------------------------------------------------------------------------------------------------------------
void 			MpdRawDataDecoder::Finish()
{
	pImpl->Finish();
}
//------------------------------------------------------------------------------------------------------------------------
void			MpdRawDataDecoder::FillEvenHeader(UInt_t runId, UInt_t eventId, UInt_t rackId)
{
	auto pFairHeader =  m_Run->GetEventHeader();
	pFairHeader->SetRunId(runId);

	auto pMpdHeader = dynamic_cast<MpdEventHeader*>(pFairHeader);
	if(pMpdHeader)
	{
		pMpdHeader->m_EventId = eventId;
		pMpdHeader->m_RackId = rackId;
	}
	else LOG(error)<<"MpdEventHeader class instance was not created. Used default FairEventHeader class."; 
}
//------------------------------------------------------------------------------------------------------------------------
bool			MpdRawDataDecoder::AddFile(const char* flnm)
{
	return pImpl->AddFile(flnm);
}
//------------------------------------------------------------------------------------------------------------------------
void 			MpdRawDataDecoder::DumpArray(const TClonesArray *array) const
{
	for(int i = 0, size = array->GetEntries(); i < size; i++) array->At(i)->Dump();
}
//------------------------------------------------------------------------------------------------------------------------

