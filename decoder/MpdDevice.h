//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_Device_H
#define __MPD_Device_H

#include <memory>
#include <iostream>
#include <list>
#include <Rtypes.h>
#include <TNamed.h>

#include "MpdRawDataDecoder.h"
//------------------------------------------------------------------------------------------------------------------------
struct MpdDeviceChannel 
{
        UInt_t sectorId = 0; 
	UInt_t detectorId = 0;
	UInt_t stripId = 0;
	UInt_t sideId = 0;

	MpdDeviceChannel(UInt_t sec, UInt_t det, UInt_t strip, UInt_t sd) : sectorId(sec), detectorId(det), stripId(strip), sideId(sd) { }
	MpdDeviceChannel() = default;

	void 		Dump(const char* comment = nullptr, std::ostream& os = std::cout) const;
};

class MpdDeviceParameters 
{
	Double_t		*m_INLs = nullptr;
	MpdDeviceChannel	*m_Channels = nullptr;	
	size_t 			nChannels = 0;
	size_t 			nBins = 0;	
public:
	MpdDeviceParameters(size_t channels, size_t bins);
	~MpdDeviceParameters();

	Double_t		GetINL(size_t channel, size_t bin) const;
	void			SetINL(size_t channel, size_t bin, Double_t value);

	MpdDeviceChannel	GetMapping(size_t channel)const;
	void			SetMapping(size_t channel, const MpdDeviceChannel&);	

	void 			DumpMapping(const char* comment = nullptr, std::ostream& os = std::cout) const;
	void 			DumpINLs(const char* comment = nullptr, std::ostream& os = std::cout) const;
};

class TClonesArray;
//------------------------------------------------------------------------------------------------------------------------
class IMpdDevice : public TNamed // interface (pure abstract class)
{
	static std::list<std::shared_ptr<IMpdDevice> >	m_list; 	// list of active devices

public:
	static MpdRawDataDecoder::RAW_DATA_FORMATS	m_dataFormat;	// default = TYPE_2022

protected:
	constexpr static size_t 	kWORDSIZE = sizeof (UInt_t); 	// 4 = 32 bits
	UInt_t				m_verboseLevel = 1;
	UInt_t				m_deviceId = (UInt_t) -1; 	// tag from TLV
	Double_t			m_slope = 0.;
	std::shared_ptr<TClonesArray>	m_container;			// data container
	size_t				m_nmbBlocks = 0;		// number of processed blocks
public:
    	IMpdDevice(std::shared_ptr<TClonesArray> array, UInt_t verbose) : m_container(array), m_verboseLevel(verbose){};
    	virtual ~IMpdDevice(){};
	virtual bool		ProcessDeviceBlock(const UInt_t *base, size_t length, UInt_t eventId, UInt_t serial) = 0;
	virtual void 		Dump(const char* comment = nullptr, std::ostream& os = std::cout) const;

	static	void		AddDevice(std::shared_ptr<IMpdDevice> device) { m_list.push_back(device);} ;
	static  bool		ProcessDeviceBlock(UInt_t m_deviceId, const UInt_t *base, size_t length, UInt_t eventId, UInt_t serial, UInt_t verbose = 1);
 	static	void 		Print(const char* comment = nullptr, std::ostream& os = std::cout);

	UInt_t			GetDeviceId()const { return m_deviceId; };
	Double_t		GetSlope()const { return m_slope; };
};
#endif 
//------------------------------------------------------------------------------------------------------------------------

