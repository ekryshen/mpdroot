//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_RAWDATA_DECODER_H
#define __MPD_RAWDATA_DECODER_H 1
//------------------------------------------------------------------------------------------------------------------------
/// \class MpdRawDataDecoder
/// 
/// \brief 
/// \author Victor Baryshnikov
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------
class FairRunAna;
#include "FairTask.h"
//------------------------------------------------------------------------------------------------------------------------
class MpdRawDataDecoder : public FairTask 
{
public:
        class Impl;
typedef	enum RAW_DATA_FORMATS 
{ 
	TYPE_OLD, 
	TYPE_2022 
} RAW_DATA_FORMATS;

typedef enum SAVE_BRANCH_FLAGS
{
        SAVE_TDC_DIGITS		= 0x1,
        SAVE_TTVXS_DIGITS	= 0x2,
        SAVE_SYNC_DIGITS	= 0x4,
        SAVE_ALL_DIGITS		= (0x1 | 0x2 | 0x4)
} SAVE_BRANCH_FLAGS;

private:
        std::unique_ptr<Impl> 		pImpl;

	std::shared_ptr<TClonesArray>	aSync;
	std::shared_ptr<TClonesArray>	aTdc;
	std::shared_ptr<TClonesArray>	aTtvxs;
	
	SAVE_BRANCH_FLAGS		fSaveFlags;
	FILE 				*fDataRawFile; // raw.dat file with raw data from DAQ
	FairRunAna* const 		m_Run = nullptr;
public:
	MpdRawDataDecoder(FairRunAna*, const char *name = "MPD RawData Decoder", Int_t verbose = 1, SAVE_BRANCH_FLAGS flags = SAVE_ALL_DIGITS, RAW_DATA_FORMATS type = TYPE_2022);
	virtual ~MpdRawDataDecoder();

	virtual InitStatus	Init();
	virtual void		Exec(Option_t * option);
	virtual void		Finish();

	bool		AddFile(const char* flnm);
	void 		DumpArray(const TClonesArray*) const;
	FairRunAna*	GetRun(){ return m_Run;}	
	void		FillEvenHeader(UInt_t runId, UInt_t eventId, UInt_t rackId);

	ClassDef(MpdRawDataDecoder,2) 
};
//------------------------------------------------------------------------------------------------------------------------
#endif 

