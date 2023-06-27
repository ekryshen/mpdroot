//------------------------------------------------------------------------------------------------------------------------
/// \class MpdRawDataDecoder
/// 
/// \brief 
/// \author Victor Baryshnikov
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------
#include <iostream>
#include <algorithm>
#include <set>

#include "FairRunAna.h"

#include "UniDbDetectorParameter.h"
#include "TLVBlock.h"
#include "MpdDevice.h"
#include "MpdDecoderImpl.h"

enum eState
{
kRunClosed,
kRunOpen
};

//------------------------------------------------------------------------------------------------------------------------
class LChain
{
typedef std::list<std::pair<FILE*, std::string> >T_fileChain;
    	T_fileChain				m_InputFiles;
	T_fileChain::iterator			m_currentFile = m_InputFiles.end();
	size_t					m_position = 0;				// [byte] base position at current open input file
	eState					m_state = kRunClosed;

public:
	//-------------------------------------------------------------------------------------------------
	bool			AddFile(const char* flnm)
	{
		auto file = fopen(flnm, "rb");
		if (file == NULL) 
		{
	        	std::cerr<<"\n[ERROR] [MpdRawDataDecoder::Impl::AddFile] - Can't open file: "<<flnm;
	        	return false;
	    	}
		m_InputFiles.push_back(std::make_pair(file, std::string(flnm)));

	return true;
	}
	//-------------------------------------------------------------------------------------------------// call AFTER all AddFile() methods
	InitStatus			Init() 		
	{
		if(m_InputFiles.empty())
		{
			cerr<<"\n[ERROR] [MpdRawDataDecoder::Impl::Init] - Input raw data files chain is empty.";
			return kFATAL;	
		}

		m_currentFile = m_InputFiles.begin();
		m_position = 0;
		m_state = kRunClosed;

	return kSUCCESS;
	}
	//-------------------------------------------------------------------------------------------------// call methods AFTER Init()						
	bool 			OpenNextFile()
	{
		if(++m_currentFile != m_InputFiles.end()) // open next file
		{
			m_position = 0;
			m_state = kRunClosed;
			return true;
		}
		return false;
	}
	//-------------------------------------------------------------------------------------------------
	eReadStatus		_checkReadError(FILE* fp) const
	{
		if(feof(fp)) return kEof;
		else if(ferror(fp)) return kReadError;
	return kOk;
	}
	//-------------------------------------------------------------------------------------------------
	void			SetRunOpen() { m_state = kRunOpen; };
	eState			GetState()const{ return m_state;};
	size_t&			GetPosition(){ return m_position;};
	size_t 			AddPosition(size_t size){ return (m_position += size);};
	const char*		GetFlnm()const{ if(m_currentFile != m_InputFiles.end()) return m_currentFile->second.c_str(); else return nullptr;};
	FILE*			GetFile(){  if(m_currentFile != m_InputFiles.end()) return m_currentFile->first; else return nullptr;};
};
//------------------------------------------------------------------------------------------------------------------------
class LTLVPipe		// fifo on 3 TLV blocks
{
	TLVBlock	*block[3] = {};
	bool		IsFirstRead = true;
	bool		IsEof = false;
public:
	LTLVPipe(Int_t verbose)
	{
		for(size_t i=0;i<3;i++)	block[i] =  new TLVBlock(verbose);
		BeginNextFile();	
	}

	void	BeginNextFile()
	{		
		block[0]->m_status = block[1]->m_status = block[2]->m_status = kUnknownError;
		IsFirstRead = true;
		IsEof = false;
	}

	void	StopReadFile(){	IsEof = true;}

	void 	Dump(const char* comment = nullptr, std::ostream& os = std::cout) const
	{
		if(comment != nullptr) os<<comment;

		block[0]->Dump(" TLV1:", os);
		block[1]->Dump(", TLV2:", os);
		block[2]->Dump(", TLV3:", os);
	}

	TLVBlock* Read(FILE* fp, size_t& position)
	{
		if(IsFirstRead)
		{
			block[0]->Read(fp, position); // read to A
			block[1]->Read(fp, position); // read to B
			block[2]->Read(fp, position); // read to C

			IsFirstRead = false;
		}
		else
		{
			std::swap(block[0], block[1]); // ABC 	-> BAC
			std::swap(block[1], block[2]); // 	-> BCA

			if(!IsEof) block[2]->Read(fp, position); 
		}
		return block[0]; // return B, read to A
	}

	eReadStatus				GetStatus(size_t i)const{ assert(i < 3); return block[i]->m_status;}
	std::pair<eReadStatus, eReadStatus>	GetStatus()const{ return std::make_pair(block[0]->m_status, block[2]->m_status);}
	UInt_t					GetTag(size_t i)const { assert(i < 3); return block[i]->m_tag;}
};
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
MpdRawDataDecoder::Impl::Impl(MpdRawDataDecoder *owner, RAW_DATA_FORMATS type, Int_t verbose)
: m_Owner(owner), m_verbose(verbose), m_Format(type)
{
	m_Input = new LChain;
	m_pipe = new LTLVPipe(m_verbose);

	IMpdDevice::m_dataFormat = m_Format;

	if(type == TYPE_OLD) // switch markers to old format
	{
		kSYNC 		= 0x2A502A50; // Regular event block	[P*P*]
		kENDOFSPILL 	= 0x4A624A62;
	} 
	else if(type == TYPE_2022) 
	{
		kSYNC 		= 0x2A50D5AF; // Regular event block	[-OP*]
		kENDOFSPILL 	= 0x4A62B59D; //  []
	}
	// add new format markers here ...
}
//------------------------------------------------------------------------------------------------------------------------
MpdRawDataDecoder::Impl::~Impl()
{
	delete m_Input;
	delete m_pipe;
}
//------------------------------------------------------------------------------------------------------------------------
InitStatus	MpdRawDataDecoder::Impl::Init()
{
	auto status = m_Input->Init();

	if(m_verbose > 1) std::cout<<" [MpdRawDataDecoder::Impl::Init] Input raw data file: "<<m_Input->GetFlnm()<<std::endl;

 return kSUCCESS;
}
//------------------------------------------------------------------------------------------------------------------------
void 			MpdRawDataDecoder::Impl::Finish()
{
	IMpdDevice::Print();
}
//------------------------------------------------------------------------------------------------------------------------
std::pair<eReadStatus, eReadStatus>		MpdRawDataDecoder::Impl::ReadNextEvent(size_t& eventId, size_t& position)
{
	auto retvalue = std::make_pair(kOk, kOk);

	switch(m_Input->GetState())
	{
		case kRunClosed:
			m_pipe->BeginNextFile();					// switch pipe to new file
		{
			auto status = ReadStartRunBlock();				// try to open Run block
		
			if(status == kOk) m_Input->SetRunOpen(); 			// switch state to kRunOpen
			else _checkFatal(status);
		}
		case kRunOpen:
		
			retvalue = ReadRegularEventBlock();				// try to read event block
			eventId = m_EventId;
			position = m_Input->GetPosition();

			if(retvalue.first != kOk) _checkFatal(retvalue.first);
		
			if(retvalue.second == kEndRun) 					// pipe 3rd block = kEndRun 
			{
				m_pipe->StopReadFile();					// kEndRun block loaded to pipe already
				auto IsOpen = m_Input->OpenNextFile();			// try to open next file into chain
				if(!IsOpen)	retvalue.second = kFinished;		// request to STOP run
			}

		break;
	}

//	if(retvalue.second != kOk) m_pipe->Dump("\n", cerr);

return retvalue;
}
//------------------------------------------------------------------------------------------------------------------------
std::pair<eReadStatus, eReadStatus> 		MpdRawDataDecoder::Impl::ReadRegularEventBlock()
{
	auto fp = m_Input->GetFile();

	auto block = m_pipe->Read(fp, m_Input->GetPosition());
	auto retval = m_pipe->GetStatus(); // (ABC); status[0] - current used block A status, status[2] - C block status
	if(m_pipe->GetTag(2) == kENDRUN) retval.second = kEndRun; // update status, if C block is last at file	

	UInt_t  payloadLength = block->m_length;
	if(m_Format == TYPE_OLD) payloadLength += 1; 

	// parse EventId 
	m_EventId = (*(block->m_base));
	size_t offset = 1;
	m_Input->AddPosition(kWORDSIZE);

//	if(m_verbose > 2) std::cerr<<"\n event="<<m_EventId<<", offset="<<toHex(m_Input->GetPosition())<<"("<<(size_t)(m_Input->GetPosition()/1048576.)<<")";

	// parse Device Event Blocks
	while(payloadLength > offset)
	{
		offset = m_eBlock.Parse(block->m_base, offset, m_Input->GetPosition()); // return offset to next DeviceEventBlock 
		
		auto ok = IMpdDevice::ProcessDeviceBlock(m_eBlock.deviceId, m_eBlock.base, m_eBlock.length, m_EventId, m_eBlock.serial, m_verbose);
		if(!ok) // unknown deviceId
		{
			//if(m_verbose > 1) std::cout<<"\n -W- [IMpdDevice::ProcessDevice] unknown deviceId: "<<toHex(m_eBlock.deviceId)<<", length: "<<m_eBlock.length<<std::endl;
		}
	}

	// Fill EventHeader branch
	m_Owner->FillEvenHeader(m_RunId, m_EventId, m_RackId); 
return retval;
}
//------------------------------------------------------------------------------------------------------------------------
eReadStatus		MpdRawDataDecoder::Impl::ReadStartRunBlock() // possible return values: kBadMarker, kEof, kReadError, kOk
{
	auto fp = m_Input->GetFile();

	auto block = m_pipe->Read(fp, m_Input->GetPosition());
	if(block->m_tag != kSTARTRUN) return kBadMarker;

	// parse inserted TLV blocks and fill MpdRunHeader
	UInt_t runNmb = -1;
	TString indexline = {};
	size_t length, offset = 0;
	UInt_t tag; const UInt_t *base = nullptr ;

	while(block->m_length > offset)
	{
		offset = block->ParseInnerTLV(offset, tag, length, &base, m_Input->GetPosition());

		switch(tag)
		{
			case kRUNNUMBER:
				runNmb = (*base);
			break;

			case kRUNINDEX: 
			{	
				UintTo4Char conv;
				size_t i = 0;
				while(i < length)
				{	
					//  uint32_t      bswap32(uint32_t int32);	 #include <sys/endian.h>
					//  uint32_t      htobe32(uint32_t host_32bits); #include <endian.h>	
					conv.in = (*(base + i));
					indexline += conv.ch;
					i++;
				}
			}
			break;

			default: // now unused TLV blocks
			break;
		}
	}

	m_RunId = runNmb;
	m_RackId = TString(indexline(15, 2)).Atoi();

return kOk;
}
//------------------------------------------------------------------------------------------------------------------------
IMpdDevice* 		MpdRawDataDecoder::Impl::CreateDevice(std::shared_ptr<IMpdDevice> device)
{
	IMpdDevice::AddDevice(device);
return device.get();
}
//------------------------------------------------------------------------------------------------------------------------
bool			MpdRawDataDecoder::Impl::AddFile(const char* flnm)
{
	return m_Input->AddFile(flnm);
}
//------------------------------------------------------------------------------------------------------------------------
void		MpdRawDataDecoder::Impl::_checkFatal(eReadStatus error)
{
	if(error == kReadError || error == kBadMarker || error == kEof)
	{
		// now unhandled errors
		std::cerr<<"\n FATAL - [MpdRawDataDecoder::Impl::ReadNextEvent]  unhandled error.";
		exit(-1);
	}
}
//------------------------------------------------------------------------------------------------------------------------

