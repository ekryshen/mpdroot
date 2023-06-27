//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_TOF_DIGIT_IMPL_H
#define __MPD_TOF_DIGIT_IMPL_H 1

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
#include <iomanip>

#include <TClonesArray.h>
#include <TTimeStamp.h>

//#include "MpdTofDigitProducer.h"
#include "MpdDevice.h"
#include "MpdTDCDigit.h"
#include "MpdTofDigit.h"
#include "UniDbDetectorParameter.h"

#include "MpdTofDigitProducer.h"
//------------------------------------------------------------------------------------------------------------------------
class MpdTofDigitProducer::Impl
{
typedef std::map<UInt_t, std::unique_ptr<MpdDeviceParameters> > 	Ttdcmap;

	MpdTofDigitProducer		*mOwner;
	Int_t 				m_verbose;
	Ttdcmap 			m_mTDCs;	// [key] = serial, [value] = strip mapping and INLs (MpdDeviceParameters class)
	std::map<UInt_t, TTimeStamp> 	m_mTSs; 	// TDC Time Shifts, [key] = serial
	std::multimap<UInt_t, UInt_t>	m_mmTTVXS;	// [key] = TVXS serial, [value] = TDC serial 
	bool				m_IsReady = false;
public:
	Impl(MpdTofDigitProducer *owner, Int_t verbose) : mOwner(owner), m_verbose(verbose) { }

	InitStatus	Init(const TOF_DIGIT_PRODUCER_DESC& desc);
	void 		Exec(const TOF_DIGIT_PRODUCER_DESC& desc);

	bool		IsReady() const {return m_IsReady; }
	Bool_t 		setMappingFromFile(const char* flnm,  size_t nChannels = 72, size_t nBins = 1024); // for TDC72VXS
	Bool_t 		setMappingFromDatabase(UInt_t req_period, UInt_t req_run, Int_t value_key, size_t nChannels = 72, size_t nBins = 1024); // for TDC72VXS

	Bool_t 		setINLsFromFile(const char* flnm, size_t nChannels = 72, size_t nBins = 1024); // for TDC72VXS
	Bool_t 		setINLsFromDatabase(const char* paramName, int period_number, int run_number, size_t nChannels = 72, size_t nBins = 1024); // for TDC72VXS
	Bool_t 		saveINLsToFile(const char* flnm, UInt_t TDCSerial, size_t nChannels = 72, size_t nBins = 1024); // for TDC72VXS

	void 		FillTimeShiftsFromSYNC();
	void 		FillTimeShiftsFromTTVXS();
	Bool_t 		SetMappingTTVXSFromFile(const char* flnm);
	void		NormalizeTimeShifts(); 

	void 		DumpMapping(const char* comment = nullptr, std::ostream& os = std::cout) const;
	void 		DumpINLs(const char* comment = nullptr, std::ostream& os = std::cout) const;
//------------------------------------------------------------------------------------------------------------------------
};
#endif 

