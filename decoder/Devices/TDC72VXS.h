#ifndef __MPD_TDC72VXS_H
#define __MPD_TDC72VXS_H

#include <TClass.h>

#include "TLVBlock.h"

#include "MpdTDCDigit.h"
#include "MpdSyncDigit.h"
#include "MpdDevice.h"
//------------------------------------------------------------------------------------------------------------------------
class TDC72VXS : public IMpdDevice 
{
	std::shared_ptr<TClonesArray>	m_syncs;

public:
	TDC72VXS(std::shared_ptr<TClonesArray> array = nullptr, std::shared_ptr<TClonesArray> sync = nullptr, UInt_t verbose = 1) : IMpdDevice(array, verbose)
	{ 
		m_deviceId = 0xD0; 
		m_slope = 24. / 1024.; // 24 ns = 1024 bins
		SetName("TDC72VXS");
		m_syncs = sync;
	};

    	virtual ~TDC72VXS() = default;
	//------------------------------------------------------------------------------------------------------------------------
	bool		ProcessDeviceBlock(const UInt_t *base, size_t length, UInt_t eventId, UInt_t serial)
	{
	        // https://afi.jinr.ru/TDC72VXS_DataFormat
		
		if(m_dataFormat == MpdRawDataDecoder::TYPE_2022)
		{
			// Parse MStream Block Header. https://afi.jinr.ru/MpdDeviceRawDataFormat
             		UInt_t data = (*base); 
			UInt_t subTypeBits	= (data & 0xFF000000) >> 24; 	// word# 0, 	24:31 Subtype-defined bits
			UInt_t blocklength	= (data & 0x00FFFFFC) >> 2; 	// 		 2:23 MStream Payload Length: S 32-bit words
			UInt_t subtype		= (data & 0x00000003); 		//  		 0:1  MStream Subtype

			// Parse WR time for TDC72_VXS
			// https://afi.jinr.ru/TAI/64-bit/Timestamp
	       		UInt_t ts_t0_s 	= (*(base + 1)); 			// word# 4?1, 	0:31 Event timestamp, seconds. 
	        	data 		= (*(base + 2)); 			// word# 5?2, 
    			UInt_t ts_t0_ns = (data & 0x00FFFFFC) >> 2;  		// 		2:31 Event timestamp, nanoseconds.	  
    			UInt_t TaiFlag  = (data & 0x00000003);   		// 		0:1  TAI flags.
	
    			auto digit = new((*m_syncs)[m_syncs->GetEntriesFast()]) MpdSyncDigit(serial, eventId, ts_t0_s, ts_t0_ns);
			if(m_verboseLevel > 2) digit->Dump();

			// Parse TDC72_VXS Data Block
	        	data = (*(base + 3)); 					// word# 5?3, 	TDC Data Block header. 	
	        	UInt_t payload 	= (data & 0x0000FFFF) / kWORDSIZE; 	//  		0:15 Data Payload length in words. 
			UInt_t blockBits= (data & 0xF0000000) >> 28; 		// 		16:27 Data Block specific bits
	        	UInt_t dataType = (data & 0xF0000000) >> 28; 		// 		28:31 Data Type. 0x0 - TDC Data, 0xF - Statistic data block
		
			if(m_verboseLevel > 3)
			{
				std::cout<<"\n [TDC72VXS::ProcessDevice] Data Block header="<<toHex(data)<<", payload="<<payload<<", dataType="<<dataType;
			}

			if(dataType == 0x0) // TDC Data
			{
				Int_t offset = 0;
				while (offset < payload) 
				{
					UInt_t data = *(base + 4 + offset++);
					UInt_t type = (data & 0xF0000000) >> 28;		
        				
					if(type == 0b0100 || type == 0b0101) // 4 - leading, 5 - trailing   28 bit 0 or 1 
					{
            					UInt_t time 	= (data & 0x001FFFFF);		//  0:20
           					UInt_t channel 	= (data & 0x0FE00000) >> 21; 	// 21:27

            					auto digit = new((*m_container)[m_container->GetEntriesFast()]) MpdTDCDigit(serial, m_deviceId, (type == 0b0100), channel, time);
						if(m_verboseLevel > 2) digit->Dump();
					}
					else if(type == 6) // TDC error
					{
						if(m_verboseLevel > 1)
						{
        						fprintf(stderr, "Warning: TDC (serial 0x%08X tdcID %d) type: %d, error code: 0x%04X\n", serial, ((data >> 24) & 0xF), type, (data & 0x7FFF)); 
						}

						return false;
					}
				}
			}
			else if(dataType == 0xF) // TDC Static Data block
			{
	                    //skip it now
	        	}
		}
		else if(m_dataFormat == MpdRawDataDecoder::TYPE_OLD)
		{
			std::cout<<"\n [TDC72VXS::ProcessDevice] MpdRawDataDecoder::RAW_DATA_FORMATS::TYPE_OLD format not implemented now.";
		}

		return true;
	}
	//------------------------------------------------------------------------------------------------------------------------
	void 		Dump(const char* comment = nullptr, std::ostream& os = std::cout) const
	{
   		IMpdDevice::Dump(comment, os);

		os<< ", classname="<<toColor::green()<<m_syncs->GetClass()->GetName()<<toColor::nocolor();
	}
};
#endif 
//------------------------------------------------------------------------------------------------------------------------

