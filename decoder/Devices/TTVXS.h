#ifndef __MPD_TTVXS_H
#define __MPD_TTVXS_H

#include <iostream>

#include "MpdTDCDigit.h"
#include "MpdDevice.h"
//------------------------------------------------------------------------------------------------------------------------
class TTVXS : public IMpdDevice 
{
	std::shared_ptr<TClonesArray>	m_syncs;

public:
	TTVXS(std::shared_ptr<TClonesArray> array, UInt_t verbose = 1) : IMpdDevice(array, verbose)
	{ 
		m_deviceId = 0xCF; 
		SetName("TTVXS");
	};

    	virtual ~TTVXS() = default;
	//------------------------------------------------------------------------------------------------------------------------
	bool		ProcessDeviceBlock(const UInt_t *base, size_t length, UInt_t eventId, UInt_t serial)
	{
/*	
		size_t offset = 0;
 		// M-Stream 2.2 protocol using Data Subtype 0.
		// parse M-Stream Header, https://afi.jinr.ru/MpdDeviceRawDataFormat
 		UInt_t data 	=  (*(base + offset++));  		// offset = 0
 		UInt_t deviceID 	= (data & 0xFF000000) >> 24;	// 24:31
 		UInt_t flags 		= (data & 0x00FC0000) >> 18;	// 18:23
 		UInt_t subtype 		= (data & 0x00030000) >> 16;	// 16:17
 		UInt_t fragmentLenth 	= (data & 0x0000FFFF);		//  0:15

	
		data 		= (*(base + offset++));			// offset = 1
 		UInt_t packetId 	= (data & 0xFFFF0000) >> 16;	// 16:31
 		UInt_t fragmentOffset 	= (data & 0x0000FFFF);		//  0:15  [bytes]

std::cerr<<"\n AAAAAAAAAAA ENTER  deviceID="<<toHex(deviceID)<<", flags="<<flags<<", subtype="<<subtype
<<", fragmentLenth="<<fragmentLenth<<", packetId="<<toHex(packetId, 4)<<", fragmentOffset="<<toHex(fragmentOffset, 4);
*/
	   	UInt_t d0,  d3;
    		UShort_t TriggerType, TriggerSource;
 
    		d0 =  (*base); //d[0] MStream Header. https://afi.jinr.ru/MpdDeviceRawDataFormat
    		//https://afi.jinr.ru/TTVXS_DataFormat
    		UInt_t d1 =   (*(base + 1)); 			// 0:31 Event timestamp, seconds .
		UInt_t data = (*(base + 2));
    		UInt_t d2 = 	(data & 0xFFFFFFFC) >> 2; 	// 2:31 Event timestamp, nanoseconds.
    		UInt_t flag = 	(data & 0x00000003); 		// 0:1	TAI flags. 

 		TTimeStamp ts(d1,d2);
//std::cerr<<"\n AAAAAAAAAAA   base="<<toHex(d0)<<" d1="<<toHex(d1)<<" d2="<<toHex(d2)<<" time="<<ts.AsDouble();

  		Int_t offset = 0;
		while (offset < length)
		{
			UInt_t data = (*(base + 3 + offset));
	        	UInt_t payload 	= (data & 0x0000FFFF) / kWORDSIZE; 	//  0:15 Data Payload length in words. 
	        	UInt_t type 	= (data & 0xF0000000) >> 28; 		// 28:31 Data Type. 0xA - Trigger Data Block, 0xF - Statistic Data Block
	
//std::cerr<<"\n AAAAAAAAAAA data="<<toHex(data,4)<<", payload="<<payload<<" type="<<toHex(type, 2);
       
	        	if(type == 0xA) 
			{
				data = (*(base + 4 + offset));
				TriggerType  = (data & 0x0000FF00) >> 8;
				TriggerSource = data & 0x000000FF;

            			offset += ((payload + 1)* kWORDSIZE); // skip Auxiliary counters
			} 
			else if(type == 0xF)
			{
        			//skip it
				offset += (payload*kWORDSIZE);
        		}
        		
if(payload < 2)exit(-1);//   LSP123 ???? 
		}
   
    		auto digit = new((*m_container)[m_container->GetEntriesFast()]) MpdTTVXSDigit(serial, ts, TriggerType, TriggerSource);
		if(m_verboseLevel > 2) digit->Dump();

		return true;
	}
};
#endif 
//------------------------------------------------------------------------------------------------------------------------

