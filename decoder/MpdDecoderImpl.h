//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_DECODER_IMPL_H
#define __MPD_DECODER_IMPL_H 1

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
#include <list>
#include <cstdarg>

#include <TClonesArray.h>
#include <TFile.h>

#include "MpdTDCDigit.h"
#include "TLVBlock.h"
#include "MpdDevice.h"
#include "MpdRawDataDecoder.h"
//------------------------------------------------------------------------------------------------------------------------
union UintTo4Char
{
    UInt_t in;
    Char_t ch[4];
};

class LChain;
class LTLVPipe;
//------------------------------------------------------------------------------------------------------------------------
class MpdRawDataDecoder::Impl
{
	MpdRawDataDecoder			*m_Owner = nullptr;
	LChain					*m_Input = nullptr;
	LTLVPipe				*m_pipe = nullptr;
	RAW_DATA_FORMATS			m_Format;

	UInt_t					kSYNC, kENDOFSPILL;
	UInt_t 					m_EventId = (UInt_t) -1;
	UInt_t 					m_RunId = (UInt_t) -1;
	UInt_t 					m_RackId = (UInt_t) -1;
	Int_t 					m_verbose;
	DeviceEventBlock			m_eBlock = {m_verbose};

	eReadStatus				_checkReadError(FILE* fp) const;
	std::pair<eReadStatus, eReadStatus>	ReadRegularEventBlock();
	eReadStatus				ReadStartRunBlock();

enum eMARKERS : UInt_t 		// Sync word (32-bit integer)
{
// https://afi.jinr.ru/MpdDeviceRawDataFormat
kSTARTRUN 	= 0x72617453, // Start run block 	[Star]
kENDRUN 	= 0x706F7453, // Stop run block 	[Stop]
kRUNNUMBER 	= 0x236E7552, // Run Number Record	[Run#]
kRUNINDEX 	= 0x78646E49, // Run Index Record	[Indx]

};
public:
	Impl(MpdRawDataDecoder *owner, RAW_DATA_FORMATS type, Int_t verbose);
	~Impl();

	InitStatus				Init();
	IMpdDevice* 				CreateDevice(std::shared_ptr<IMpdDevice> device);
	bool					AddFile(const char* flnm);
	void 					Finish();
	void					_checkFatal(eReadStatus error);
	std::pair<eReadStatus, eReadStatus>	ReadNextEvent(size_t& eventId, size_t& position);

//------------------------------------------------------------------------------------------------------------------------
};
#endif 

