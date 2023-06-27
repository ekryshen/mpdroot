//------------------------------------------------------------------------------------------------------------------------
#include <iostream>

#include <TClass.h>
#include <TClonesArray.h>

#include "TLVBlock.h"
#include "MpdDevice.h"
//------------------------------------------------------------------------------------------------------------------------ 
MpdDeviceParameters::MpdDeviceParameters(size_t channels, size_t bins)
: nChannels(channels), nBins(bins)
{
	auto size = nChannels*nBins;
assert(size < 1000000); //  < 1k * 1k

	m_INLs = new Double_t[size];
	memset(m_INLs, 0, sizeof(m_INLs));

	m_Channels = new MpdDeviceChannel[nChannels];	
}
//------------------------------------------------------------------------------------------------------------------------ 
MpdDeviceParameters::~MpdDeviceParameters()
{
	delete[]  m_INLs;
	delete[]  m_Channels;
}
//------------------------------------------------------------------------------------------------------------------------ 
Double_t		MpdDeviceParameters::GetINL(size_t channel, size_t bin) const
{
assert(channel < nChannels);
assert(bin < nBins);
	
	return  *(m_INLs + channel*nBins +  bin); 
}
//------------------------------------------------------------------------------------------------------------------------ 
void			MpdDeviceParameters::SetINL(size_t channel, size_t bin, Double_t value)
{
assert(channel < nChannels);
assert(bin < nBins);
	
	*(m_INLs + channel*nBins +  bin) = value;
}
//------------------------------------------------------------------------------------------------------------------------ 
MpdDeviceChannel	MpdDeviceParameters::GetMapping(size_t channel)const
{
assert(channel < nChannels);

	return  *(m_Channels + channel);
}
//------------------------------------------------------------------------------------------------------------------------ 
void			MpdDeviceParameters::SetMapping(size_t channel, const MpdDeviceChannel& params)
{
assert(channel < nChannels);

        auto entry = (m_Channels + channel);

	entry->sectorId = params.sectorId; 
	entry->detectorId = params.detectorId; 
	entry->stripId = params.stripId; 
	entry->sideId = params.sideId;
}
//------------------------------------------------------------------------------------------------------------------------ 
void 			MpdDeviceParameters::DumpMapping(const char* comment, std::ostream& os) const
{
	if(comment != nullptr) os<<comment;

	for(size_t channel = 0; channel < nChannels; channel++)
	{
		(m_Channels + channel)->Dump("\n ", os);
	}
}
//------------------------------------------------------------------------------------------------------------------------ 
void 			MpdDeviceParameters::DumpINLs(const char* comment, std::ostream& os) const
{
	if(comment != nullptr) os<<comment;

        for(size_t channel = 0; channel < nChannels; channel++) 
	{
        	os<<"\n channel="<<channel<<"\n";
            
		size_t counter = 0;
		for(size_t bin = 0; bin < nBins; bin++)
		{
                	os<<" "<<GetINL(channel, bin);
			if(counter++ % 10 == 0) os<<std::endl;
        	}
	}
}
//------------------------------------------------------------------------------------------------------------------------ 
void 			MpdDeviceChannel::Dump(const char* comment, std::ostream& os) const
{
	if(comment != nullptr) os<<comment;

	os<<" sector="<<sectorId<<", detector="<<detectorId<<", strip="<<stripId<<", side="<<sideId;
}
//------------------------------------------------------------------------------------------------------------------------ 
std::list<std::shared_ptr<IMpdDevice> >		IMpdDevice::m_list;
MpdRawDataDecoder::RAW_DATA_FORMATS		IMpdDevice::m_dataFormat = MpdRawDataDecoder::TYPE_2022;
//------------------------------------------------------------------------------------------------------------------------ 
bool		IMpdDevice::ProcessDeviceBlock(UInt_t deviceId, const UInt_t *base, size_t length, UInt_t eventId, UInt_t serial, UInt_t verbose)
{
	for(auto& entry : m_list) // cycle by devices
	{
		if(entry->m_deviceId == deviceId)
		{
			entry->ProcessDeviceBlock(base, length, eventId, serial);
			entry->m_nmbBlocks++;

			if(verbose > 2) std::cout<<"\n [IMpdDevice::ProcessDeviceBlock] deviceId="<<toHex(deviceId, 2)<<", base="<<base<<", length="<<length
					<<" ("<<entry->GetName()<<") "<<entry->GetTitle();
			return true;
		}
	}

	if(verbose > 1) std::cout<<toColor::red()<<"\n [IMpdDevice::ProcessDeviceBlock] deviceId="<<toHex(deviceId, 2)<<", base="<<base
		<<", length="<<length<<" --- NOT found."<<toColor::nocolor();

return false; // guid don't found
}
//------------------------------------------------------------------------------------------------------------------------
void 		IMpdDevice::Print(const char* comment, std::ostream& os)
{
	os<<"\n [IMpdDevice::Print] -------------------------------------------------->"; 
	if(nullptr != comment) os<<comment;
	for(auto& entry : m_list) entry->Dump("\n", os);
	os<<"\n -----------------------------------------------------------------------";
}
//------------------------------------------------------------------------------------------------------------------------
void 		IMpdDevice::Dump(const char* comment, std::ostream& os) const
{
	if (comment != nullptr) os << comment;

	os<<toColor::green()<<"("<<GetName()<<") "<<toColor::nocolor()<<GetTitle()<<", deviceId="<<toColor::green()<<toHex(m_deviceId, 2)
		<<toColor::nocolor()<<", classname="<<toColor::green()<<m_container->GetClass()->GetName()<<toColor::nocolor()<<", blocks processed="<<m_nmbBlocks;
}
//------------------------------------------------------------------------------------------------------------------------
