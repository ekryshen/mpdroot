//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_TLVBlock_H
#define __MPD_TLVBlock_H

#include <iostream>
#include <iomanip>
#include <cassert>

#include <Rtypes.h>

#define assertm(exp, msg) assert(((void)msg, exp))
//------------------------------------------------------------------------------------------------------------------------
class toHex
{
	int 	value;
	int	width;
public:
	toHex(int v, int w = 8): value(v), width(w){}

	friend std::ostream& operator<<(std::ostream& os, const toHex& dt)
	{
		std::ios state(nullptr);
		state.copyfmt(os); // save current formatting

    		os<<std::hex<<std::uppercase<<"0x";
		if(dt.width>0) os<<std::setw(dt.width)<< std::setfill('0');
		os<<dt.value;

		os.copyfmt(state); // restore previous formatting
    		return os;
	}
};
//#define toHex
//------------------------------------------------------------------------------------------------------------------------
class toColor
{
	static const char* array(int index) 
	{
    		static const char* a[] = { "\033[0m", "\033[;32m", "\033[1;31m", "\033[33;44m"};
    		return a[index];
  	} 
public:
	static const char* nocolor(){ return array(0); }
	static const char* green(){ return array(1); }
	static const char* red(){ return array(2); }
	static const char* yellow(){ return array(3); }

};
//------------------------------------------------------------------------------------------------------------------------
constexpr size_t kWORDSIZE = sizeof (UInt_t); //  = 4,  DWORD = 32 bits

enum eReadStatus
{
kOk = 0,
kEof,			// feof()
kReadError,		// ferror()
kBadMarker,		// unexpected marker
kEndRun,		// current file stop run block reached
kBigBufferSize,		// requed buffer size exceeds the limit
kFinished,		// tried to read eventNmb more than exist (all files have been processed)
kUnknownError
};

//------------------------------------------------------------------------------------------------------------------------
class TLVBlock
{
	Int_t 					m_verbose;
	UInt_t 					m_2words[2]; 		// Marker + length ( 32 + 32 bytes)
	static const size_t			m_bufferSize = 10000000;
	UInt_t 					m_bigBuffer[m_bufferSize];	// big buffer for fread()
public:
	eReadStatus	m_status;
	UInt_t		m_tag;
	size_t		m_length;  		// [] = UInt_t  (*4 bytes)
	const UInt_t	*m_base = m_bigBuffer;	// first UInt_t of data

    	TLVBlock(Int_t verbose): m_verbose(verbose){};
    	virtual ~TLVBlock() = default;

	//------------------------------------------------------------------------------------------------------------------------
	// [offset] = UInt_t  (*4 bytes), return offset to next TLV block
 	size_t 		ParseInnerTLV(size_t offset, UInt_t& tag, size_t& length, const UInt_t **base, size_t& positon) const 
	{
		tag = *(m_base + offset);
		length = (*(m_base + offset + 1))/ kWORDSIZE;
		(*base) = (m_base + offset + 2);

		if(m_verbose > 1)
		{
			std::cout<<"\n [TLVBlock::ParseInnerTLV] tag="<<toColor::green()<<toHex(tag)<<toColor::nocolor()<<", position="<<toHex(positon,-1);
		}

		positon += ((2+ length) *kWORDSIZE);

	return offset + 2 + length;
	}
	//------------------------------------------------------------------------------------------------------------------------
	eReadStatus	Read(FILE* fp, size_t& positon)
	{

///std::cerr<<"\n 11"<<" feof="<<feof(fp)<<" ferr="<<ferror(fp);
		auto size = fread(m_2words, kWORDSIZE, 2, fp); 
		if(size != 2) return _checkReadError(fp, kUnknownError);

///std::cerr<<"\n 22 size="<<size<<" feof="<<feof(fp)<<" ferr="<<ferror(fp);

		m_tag = m_2words[0];
		m_length = m_2words[1] / kWORDSIZE;


		if(m_length >= 100000)// what the constant? ???????????????????????????????????????????????????? LSP
		{ 
			std::cerr<<"\n [TLVBlock::Read] Wrong data event block size: "<<m_length<<",  skip this event.\n";
			return (m_status = kBigBufferSize);
        	}

//		if(m_verbose > 1)
//		{
//			auto color = (expectedTag == m_tag) ? toColor::green() : toColor::red();
//			std::cout<<"\n [TLVBlock::Read] tag="<<color<<toHex(m_tag)<<"("<<tagName<<")"<<toColor::nocolor()
//					<<" payloadLength="<<m_length<<"("<<toHex(m_2words[1])<<"), position="<<toHex(positon,-1);
//		}

		positon += (2*kWORDSIZE);

assertm(m_length < m_bufferSize, " TLB block length > inner buffer size.");

///std::cerr<<"\n 33"<<" feof="<<feof(fp)<<" ferr="<<ferror(fp);
		// read event block to buffer	
		size = fread(m_bigBuffer, kWORDSIZE, m_length, fp);
		if(size != m_length) return _checkReadError(fp, kUnknownError);	
	
///std::cerr<<"\n 44 size="<<size<<" feof="<<feof(fp)<<" ferr="<<ferror(fp);

	return (m_status = kOk);
	}
	//------------------------------------------------------------------------------------------------------------------------
	eReadStatus		_checkReadError(FILE* fp, eReadStatus status) 
	{
		if(feof(fp)) return (m_status = kEof);
		else if(ferror(fp)) return (m_status = kReadError);
	return status;
	}
	//------------------------------------------------------------------------------------------------------------------------
	void 		Dump(const char* comment = nullptr, std::ostream& os = std::cout) const
	{
		if(comment != nullptr) os<<comment;
		os<<"(T="<<toHex(m_tag)<<", S="<<m_status<<", L="<<m_length<<")";
	}
};
//------------------------------------------------------------------------------------------------------------------------
struct DeviceEventBlock
{
	Int_t 	m_verbose;
	UInt_t	serial;
	UInt_t	deviceId;
	size_t	length;  			// [] = UInt_t  (*4 bytes)
	const UInt_t	*base = nullptr;	// payload base address

	DeviceEventBlock(Int_t verbose): m_verbose(verbose){};

	size_t 		Parse(const UInt_t *base, size_t offset, size_t& position) // [offset] = UInt_t  (*4 bytes)
	{
		serial = *(base + offset);	 	 	// 0:31
		deviceId = (*(base + offset + 1)) >> 24; 	// 24:31
		length = (*(base + offset + 1)) & 0x00FFFFFF; 	// 0:23
		length /= kWORDSIZE; // [byte] -> [DWORD]
		this->base = (base + offset + 2);		

		if(m_verbose > 1)
		{
			std::cout<<"\n [DeviceEventBlock::Parse] serial="<<toColor::yellow()<<toHex(serial)<<toColor::nocolor()
					<<", deviceId="<<toColor::yellow()<<toHex(deviceId, 2)<<toColor::nocolor()<<", length="<<length
					<<"("<<toHex(*(base + offset + 1))<<"), offset="<<offset<<", position="<<toHex(position, -1);
		}

		position += ((2+ length) *kWORDSIZE); // [DWORD] -> [byte] 

	return offset + 2 + length;
	}
};

#endif 
//------------------------------------------------------------------------------------------------------------------------

