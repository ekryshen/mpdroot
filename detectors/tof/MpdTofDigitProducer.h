//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_TOF_DIGIT_PRODUCER_H
#define __MPD_TOF_DIGIT_PRODUCER_H 1

//------------------------------------------------------------------------------------------------------------------------
/// \class MpdRawDataDecoder
/// 
/// \brief 
/// \author Victor Baryshnikov
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------

#include <iostream>
#include <TObject.h>
#include <string>

#include <TFile.h>

#include "MpdTofDigit.h"
#include "MpdSyncDigit.h"

#include "FairTask.h"

typedef enum TOF_DIGIT_TIMESHIFT
{
	TS_TTVXS,
	TS_SYNC,
	TS_NONE

} TOF_DIGIT_TIMESHIFT;

typedef struct TOF_DIGIT_PRODUCER_DESC
{
	TOF_DIGIT_TIMESHIFT	TS_mode = TS_TTVXS;	//	Time shift source = TTVXS (default)
	std::string		INL_flnm = {};		//	INL loaded from file(1), if flnm not empty; otherwise from dB(2) (default)	
	std::string		stripMap_flnm = {};	//	Strip <-> TDC channel map loaded from file(3), if flnm not empty; otherwise from dB(4) (default)
	std::string		TTVXSMap_flnm = {};	//	TTVXS serial -> TDC serial mapping
	Int_t 			periodId = 0;		// 	dB request parameter, used by method (4)
	Int_t 			runId = -1;		// 	dB request parameter, used by method (4)
	Int_t 			rackId = -1;		// 	dB request parameter, used by method (4)
}	TOFDIGIT_PRODUCER_DESC;
//------------------------------------------------------------------------------------------------------------------------
class MpdEventHeader;
class FairRunAna;
class MpdTofDigitProducer: public FairTask 
{
public:
        class Impl;

	TClonesArray 		*aTofDigits = nullptr; 

        TClonesArray 		*aSync = nullptr; 
        TClonesArray 		*aTtvxs = nullptr; 
        TClonesArray 		*aTdc = nullptr; 

private:
        std::unique_ptr<Impl> 		pImpl;

	TOF_DIGIT_PRODUCER_DESC	m_desc;
	std::string		m_flnm = {};
	TTree 			*m_Tree = nullptr;
	FairRunAna* const 		m_Run = nullptr;

public:
        MpdTofDigitProducer();
	MpdTofDigitProducer(FairRunAna*, TOF_DIGIT_PRODUCER_DESC& desc, Int_t verbose = 0, const char *name = "TofDigitProducer");
       
	virtual ~MpdTofDigitProducer();

	virtual InitStatus	Init();
	virtual void		Exec(Option_t * option);
        virtual void		Finish();
		
	// Add input file with device digits(MpdTDCDigit, MpdSyncDigit, MpdTTVXSDigit etc.). At this case, MUST BE FIRST FairTask at FairRunAna chain.
	void	AddInputFile(const char* flnm) { m_flnm = flnm; };	

ClassDef(MpdTofDigitProducer,1) 
};
//------------------------------------------------------------------------------------------------------------------------
#endif 

