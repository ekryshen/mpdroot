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

#include "MpdTDCDigit.h"
#include "MpdSyncDigit.h"
#include "MpdTTVXSDigit.h"

#include "MpdTofDigit.h"
#include "UniDbDetectorParameter.h"

#include "TDC72VXS.h"

#include "MpdTofDigitImpl.h"

// Comparator of TDCDigits by strip guid (MpdTDCDigit::operator<) 
struct lessPtrByGuid
{
	bool operator()(const MpdTDCDigit*  lhs, const MpdTDCDigit*  rhs) const 
	{
		return (*lhs) < (*rhs); 
	}
};
//------------------------------------------------------------------------------------------------------------------------
InitStatus	MpdTofDigitProducer::Impl::Init(const TOF_DIGIT_PRODUCER_DESC& desc)//, Int_t runNmb, Int_t periodId)
{
	if(m_IsReady) return kSUCCESS;

	// Mapping strip <-> TDC channel 
	bool status;
	if(desc.stripMap_flnm.empty()) // loaded from dB(4)
	{

cerr<<"\n setMappingFromDatabase  period="<<desc.periodId<<", run="<<desc.runId<<", value_key="<<desc.rackId;

		status = setMappingFromDatabase(desc.periodId, desc.runId, desc.rackId); 
	}
	else				// loaded from file(3)
	{
		status = setMappingFromFile(desc.stripMap_flnm.c_str()); 
	}

        if(status == false)
	{
            cerr<<"\n [MpdTofDigitProducer::Impl::Init] Error loading of the strip to TDC channel mapping.\n";
            return kERROR;
        }

	// Load INL parameters for mapped TDC
	if(desc.INL_flnm.empty()) 	// loaded from dB(2)
	{
		status = setINLsFromDatabase("TOF1_inl", desc.periodId, desc.runId); 
	}
	else				// loaded from file(1)
	{
		status = setINLsFromFile(desc.INL_flnm.c_str());
	}

        if(status == false)
	{
            cerr<<"\n [MpdTofDigitProducer::Impl::Init] Error loading of INL parameters.\n";
            return kERROR;
        }

	// Load Time Shift mapping for TTVXS mode
	if(desc.TS_mode == TS_TTVXS)
	{
		status = SetMappingTTVXSFromFile(desc.TTVXSMap_flnm.c_str());
	
	        if(status == false)
		{
	            cerr<<"\n [MpdTofDigitProducer::Impl::Init] Error loading Time Shift mapping for TTVXS mode.\n";
	            return kERROR;
	        }
	}

 	LOG(info)<<"[MpdTofDigitProducer::Init] (periodId="<<desc.periodId<<", runId="<<desc.runId<<")Initialization finished succesfully.";

	m_IsReady = true;

return kSUCCESS; 
}
//------------------------------------------------------------------------------------------------------------------------
void 		MpdTofDigitProducer::Impl::Exec(const TOF_DIGIT_PRODUCER_DESC& desc)
{
static double slopeTDC72VXS = TDC72VXS{}.GetSlope();
static const  UInt_t deviceId = TDC72VXS{}.GetDeviceId();

	// Fill TDC time shifts
	if(	desc.TS_mode == TS_TTVXS) 	FillTimeShiftsFromTTVXS();
	else if(desc.TS_mode == TS_SYNC) 	FillTimeShiftsFromSYNC();

	// Sort MpdTDCDigit pointers by suid
	multiset<MpdTDCDigit*, lessPtrByGuid> sortedDigits; 
	for(size_t i = 0, size = mOwner->aTdc->GetEntriesFast(); i < size; i++)
	{
		auto digit = (MpdTDCDigit*) mOwner->aTdc->At(i);
		
		if(digit->GetType() == deviceId)
		{
			sortedDigits.insert(digit);
		}
	}


//cerr<<"\n AAA aTdc size="<<mOwner->aTdc->GetEntries()<<" "<<mOwner->aTdc->GetEntriesFast();;
size_t i0=0, i1 = 0, i2 = 0, i3 =0;

	for(auto it1 = sortedDigits.cbegin(), itEnd = sortedDigits.cend(); it1 != itEnd; it1 = sortedDigits.upper_bound(*it1)) // cycle by same suid MpdTDCDigit
	{

		size_t counter = sortedDigits.count(*it1);

		if(counter < 2) continue;	// must be >= 2 MpdTDCDigit for same suid

		UInt_t serial = (*it1)->GetSerial();
		UShort_t channel = (*it1)->GetChannel();

		UInt_t value1 = 0, value2 = 0;
		Int_t leading1 = -1, leading2 = -1;

		if(counter == 2) // ideal case	
		{
			auto it2 = std::next(it1, 1);
assert(serial == (*it2)->GetSerial());
assert(channel == (*it2)->GetChannel());

			leading1 = (*it1)->GetLeading();		
			leading2 = (*it2)->GetLeading();
			value1 = (*it1)->GetValue();
			value2 = (*it2)->GetValue();
		}
		else
		{

//cerr<<"\n AA  ----------------------counter="<<counter;

			// Find MpdTDCDigit pair(leading, tailing) with minimum fValue
			auto range = sortedDigits.equal_range(*it1);
			for(auto it = range.first; it != range.second; it++)
			{
				UInt_t value = (*it)->GetValue();
	
				if((*it)->GetLeading())	// leading digit
				{	
					leading1 = 1;			
					if(value1 == 0)  	value1 = value; // first time
					else if(value < value1)	value1 = value; // new fastest found
				}
				else			// tailing digit
				{
					leading2 = 0;			
					if(value2 == 0)  	value2 = value; // first time
					else if(value < value1)	value2 = value; // new fastest found
				}



//cerr<<"\n AAA  value="<<value<<"  leading="<<(*it)->GetLeading()<<"  1="<<value1<<" 2="<<value2;

			}
		}
		
		if(leading1 + leading2 != 1) // must be 0+1 or 1+0
		{
			if(m_verbose > 3) cout<<"\n[MpdTofDigitProducer::Impl::Exec] invalid digit types: type1="<<leading1<<", type2="<<leading2;
//cerr<<"\nAAAAA invalid digit types: type1="<<leading1<<", type2="<<leading2;
			continue;
		}
	
		// Get INLs	
		Double_t INL1 = 0., INL2 = 0.;
		auto iter = m_mTDCs.find(serial);
		if(iter != m_mTDCs.end())
		{
		 	INL1 = iter->second->GetINL(channel, value1 % 1024);
		 	INL2 = iter->second->GetINL(channel, value2 % 1024); 

/////cerr<<"\n INL1="<<INL1<<"("<<value1<<", "<<leading1<<"), INL2="<<INL2<<"("<<value2<<", "<<leading2<<")";

		}
		else
		{
//cerr<<"\n TDC serial="<<hex<<serial<<dec<<" INLs don't found.";

			if(m_verbose > 0) cout<<" [MpdTofDigitProducer::Impl::Exec] TDC serial="<<hex<<serial<<dec<<" INLs don't found.\n";
			continue; // 12345LSP REMOVE!!!!
		}
	
		Double_t time1 = (value1 + INL1) * slopeTDC72VXS; 
		Double_t time2 = (value2 + INL2) * slopeTDC72VXS; 

///cerr<<"\n time1="<<time1<<", time2="<<time2<<" slope="<<slopeTDC72VXS;
		
		// Select time1 to leading time
		if(! leading1) std::swap(time1, time2);

		if(time1 > time2)
		{
/////////12345LSP			cerr<<"[MpdTofDigitProducer::Impl::Exec] Error: negative dt, time1="<<time1<<"("<<value1<<") > time2="<<time2<<"("<<value2<<").\n";
			continue;
		}

		// Get TDC time shift
		Double_t timeShift = 0.;
		if(desc.TS_mode != TS_NONE)
		{
			auto it = m_mTSs.find(serial);
			if(it == m_mTSs.end())
			{
				cerr<<"\n[MpdTofDigitProducer::Impl::Exec] Error: TimeShift for TDC serial="<<hex<<serial<<dec<<" don't found.";
				continue;
			}

			if(it->second.GetSec() == 0) timeShift = it->second.GetNanoSec();
            		else cerr<<"\n[MpdTofDigitProducer::Impl::Exec] Error: TimeShift ~ sec ("<<it->second.GetSec()<<").";
		}

		// Get mapping
		const auto& entry = iter->second->GetMapping(channel);

		// Add MpdTofDigit
		auto digit = new((*mOwner->aTofDigits)[mOwner->aTofDigits->GetEntries()]) 
				MpdTofDigit(entry.sectorId, entry.detectorId, entry.stripId, entry.sideId, time1 + timeShift, time2 - time1);

//		if(m_verbose > 3) digit->print("\n add MpdTofDigit: ");	
	}


//cerr<<"\n AAA aTofDigits size="<<mOwner->aTofDigits->GetEntries()<<" "<<mOwner->aTofDigits->GetEntriesFast();

}
//------------------------------------------------------------------------------------------------------------------------
void 		MpdTofDigitProducer::Impl::FillTimeShiftsFromSYNC()
{
	// Cleanup prev. event time shifts
	m_mTSs.clear();

    	for(Int_t i = 0; i < mOwner->aSync->GetEntriesFast(); i++) // cycle by MpdSyncDigit
	{
        	auto digit = (MpdSyncDigit*) mOwner->aSync->At(i);

		// insert new entry
		auto retval = m_mTSs.insert(std::make_pair(digit->GetSerial(), TTimeStamp((time_t)digit->GetTime_sec(), digit->GetTime_ns())));
		if(retval.second == false)
		{
			if(m_verbose > 2) cout<<"\n[MpdTofDigitProducer::Impl::FillTimeShiftsFromSYNC] TimeShift for key(TDC serial)="<<hex<<digit->GetSerial()<<dec<<" already exist!";
		}
	}

	NormalizeTimeShifts();
}
//------------------------------------------------------------------------------------------------------------------------
void 		MpdTofDigitProducer::Impl::FillTimeShiftsFromTTVXS()
{
	// Cleanup prev. event time shifts
	m_mTSs.clear();

	// Check mapping is loaded
	if(m_mmTTVXS.empty())
	{
		if(m_verbose > 2) cout<<"\n[MpdTofDigitProducer::Impl::FillTimeShiftsFromTTVXS] TTVXS serial -> TDC serial mapping don't loaded.";
		return;
	}

    	for(Int_t i = 0; i < mOwner->aTtvxs->GetEntriesFast(); i++) // cycle by MpdTTVXSDigit
    	{
        	auto digit = (MpdTTVXSDigit*) mOwner->aTtvxs->At(i);

		auto range = m_mmTTVXS.equal_range(digit->GetSerial());
		for(auto it = range.first; it != range.second; it++) // cycle by same TTVXS serial entries
		{
			// insert new entry
			auto retval = m_mTSs.insert(std::make_pair(it->second, digit->GetTime()));
			if(retval.second == false)
			{
				if(m_verbose > 2) cout<<"\n[MpdTofDigitProducer::Impl::FillTimeShiftsFromTTVXS] TimeShift for key(TDC serial)="
							<<hex<<it->second<<dec<<" already exist!";
			}
		}
	}

	NormalizeTimeShifts();
}
//------------------------------------------------------------------------------------------------------------------------
void 		MpdTofDigitProducer::Impl::NormalizeTimeShifts()
{
	if(m_mTSs.empty())return;

	// Calc. minimum TimeShift
	TTimeStamp tmin = (*std::min_element(m_mTSs.begin(), m_mTSs.end(), [](auto a, auto b) { return a.second < b.second; })).second;

	// TTimeStamp have only Add method, invert tmin
	tmin.SetSec(-1.*tmin.GetSec());
	tmin.SetNanoSec(-1.*tmin.GetNanoSec());

	// Subtract tmin
	for(auto& entry : m_mTSs)
	{
		entry.second.Add(tmin);
	}
}
//------------------------------------------------------------------------------------------------------------------------
Bool_t 		MpdTofDigitProducer::Impl::SetMappingTTVXSFromFile(const char* flnm)
{
	ifstream fs(flnm);
	if(!fs) 
	{
        	cerr<<"\n[MpdTofDigitProducer::Impl::SetMappingTTVXSFromFile] Can't open the file"<<flnm; 
        	return false;
    	}
	 
    	while(!fs.eof()) 
    	{
		UInt_t TTVXS_serial, fcon, TDC_serial;
        	fs>>std::hex>>TTVXS_serial>>std::dec>>fcon>>std::hex>>TDC_serial>>std::dec;

        	if(m_verbose > 2) cout<<"\n serials: TTVXS="<<hex<<TTVXS_serial<<", TDC="<<TDC_serial<<dec;

		m_mmTTVXS.insert(std::make_pair(TTVXS_serial, TDC_serial));
    	}

    	fs.close();
return true;
}
//------------------------------------------------------------------------------------------------------------------------
Bool_t 		MpdTofDigitProducer::Impl::setMappingFromFile(const char* flnm, size_t nChannels, size_t nBins) 
{
	fstream file;
	file.open(flnm, std::fstream::in);

	if(file.fail()) 
	{
		cerr<<"\n [MpdTofDigitProducer::Impl::setMapFromFile] Cannot open the file \""<<flnm<<"\""; 
		return false;
	}
        else cout<<"\n [MpdTofDigitProducer::Impl::setMapFromFile] Reading TOF StripMap file \""<<flnm<<"\"";

	// cached values
	UInt_t lastkey = 0;
	Ttdcmap::iterator iter = m_mTDCs.end();

	UInt_t serial, channel, sector, detector, strip;
	char side_c;

	while(!file.eof()) 
	{
		file>>std::hex>>serial>>std::dec>>channel>>sector>>detector>>strip>>side_c;
		if(file.eof()) { break; } 

		// Fast find TDC by serial
		if(serial != lastkey) // need update cached values
		{
			iter = m_mTDCs.find(serial);

			if(iter == m_mTDCs.end()) // need create new map entry 
			{			
				iter = m_mTDCs.insert(std::make_pair(serial, std::make_unique<MpdDeviceParameters>(nChannels, nBins))).first;
			}
		
			lastkey = serial;
		}

		// Fill map value
		iter->second->SetMapping(channel, MpdDeviceChannel(sector, detector, strip, (side_c == 'L') ? 0 : 1));
	}

	file.close();

//	if(m_verbose > 4) DumpMapping("\n [MpdTofDigitProducer::Impl::setMapFromFile] loaded mapping:");

return true;
}
//------------------------------------------------------------------------------------------------------------------------
Bool_t 		MpdTofDigitProducer::Impl::setMappingFromDatabase(UInt_t req_period, UInt_t req_run, Int_t value_key, size_t nChannels, size_t nBins) 
{

req_period = 0;
//value_key = 5;

//if(m_verbose > 1) 
cerr<<"\n[MpdTofDigitProducer::Impl::setMapFromDatabase] looking for StripMap in database for period="<<req_period<<", run="<<req_run<<", value_key="<<value_key;

	UniDbDetectorParameter* dbDetParam = UniDbDetectorParameter::GetDetectorParameter("TOF1", "TOF1_StripMap", req_period, req_run, value_key);
	if (dbDetParam == nullptr) 
	{
		cerr<<"\n[MpdTofDigitProducer::Impl::setMappingFromDatabase]  Failed to load strip mapping from database!";
		return false;
	}

    	vector<UniValue*> array;
    	dbDetParam->GetValue(array);

	// dbDetParam::TOF1_StripMap structure:
	//	int serial (hex);
	//	int channel;
	//	value [0]: sector
	//	value [1]: detector
	//	value [2]: strip
	//	value [3]: side (0 = left, 1 = right)
   	cout << "Found stripmap in db, size = " << array.size() << endl;

	// cached values
	UInt_t lastkey = 0;
	Ttdcmap::iterator iter = m_mTDCs.end();

	for (Int_t i = 0; i < array.size(); i++) 
	{
		// Parse dB entry
		auto entry = (MapDVectorValue*) array.at(i); // taking single row (serial) from map
		UInt_t serial = entry->serial; 
		UInt_t channel = entry->channel;
		UInt_t sector = entry->value[0];
		UInt_t detector = entry->value[1];
		UInt_t strip = entry->value[2];
		UInt_t side = entry->value[3];	

		// Fast find TDC by serial
		if(serial != lastkey) // need update cached values
		{
			iter = m_mTDCs.find(serial);

			if(iter == m_mTDCs.end()) // need create new map entry 
			{			
				iter = m_mTDCs.insert(std::make_pair(serial, std::make_unique<MpdDeviceParameters>(nChannels, nBins))).first;
			}
		
			lastkey = serial;
		}

		// Fill map value
		iter->second->SetMapping(channel, MpdDeviceChannel(sector, detector, strip, side));       
	}

//	if(m_verbose > -1) DumpMapping("\n [MpdTofDigitProducer::Impl::setMapFromDatabase] loaded mapping:");

return true;
}
//------------------------------------------------------------------------------------------------------------------------
Bool_t 		MpdTofDigitProducer::Impl::saveINLsToFile(const char* flnm, UInt_t TDCSerial, size_t nChannels, size_t nBins)
{
	auto it = m_mTDCs.find(TDCSerial);
	if(it == m_mTDCs.end()) 
	{
		cerr<<"\n[MpdRawDataDecoder::Impl::saveINLToFile] - TDC serial="<<hex<<TDCSerial<<dec<<" isn't in the placement map.";
		return false;
	}

	fstream fs(flnm, fstream::out);

	fs<<"[TDC-"<<hex<<TDCSerial<<dec<<"-inl_corr]"<<endl;

	for(size_t channel = 0; channel < nChannels; channel++) 
	{
		fs<<channel<<"=";

		for(size_t bin = 0; bin < nBins; bin++) 
		{
			fs<<it->second->GetINL(channel, bin);
			if(bin != nBins-1) fs<<", ";		
		}

		if(channel != nChannels - 1) fs<<endl;
	}

	fs.close();

return true;
}
//------------------------------------------------------------------------------------------------------------------------
Bool_t 			MpdTofDigitProducer::Impl::setINLsFromDatabase(const char* paramName, int period_number, int run_number, size_t nChannels, size_t nBins) 
{

period_number = 0;


//if(m_verbose > 1) 
cerr<<"\n[MpdTofDigitProducer::Impl::setINLsFromDatabase] looking for INLs in database for period="
				<<period_number<<", run="<<run_number<<", paramName="<<paramName;



    	for(auto it = m_mTDCs.begin(), itEnd = m_mTDCs.end(); it != itEnd; it++)
	{
        	UInt_t TDC_serial = it->first;
        	UniDbDetectorParameter* dbDetParam = UniDbDetectorParameter::GetDetectorParameter("TOF1", paramName, period_number, run_number, TDC_serial);
        
		if(dbDetParam == nullptr) 
		{
	            	cerr<<"\n[MpdTofDigitProducer::Impl::setINLsFromDatabase] - Failed to load INLs from database! TDC_serial= "<<hex<<TDC_serial<<dec<<", param. name= "<<paramName;
            		continue;
        	}

        	vector<UniValue*> array; 
        	dbDetParam->GetValue(array);

        	auto entry = (TdcInlValue*) array.at(0);
        	for(size_t channel = 0; channel < nChannels; channel++) 
		{           
			for(size_t bin = 0; bin < nBins; bin++)
			{
                		it->second->SetINL(channel, bin, entry->inl[channel][bin]);
			}
        	}

        	if(m_verbose > 2) cout<<"\n[MpdTofDigitProducer::Impl::setINLsFromDatabase] INLs for TDC serial="<<hex<<TDC_serial<<dec<<" loaded succesfully from DB.";
	}

//	if(m_verbose > 5) DumpINLs("\n [MpdTofDigitProducer::Impl::setINLsFromDatabase] loaded INLs:");

return true;
}
//------------------------------------------------------------------------------------------------------------------------
Bool_t 			MpdTofDigitProducer::Impl::setINLsFromFile(const char* flnm, size_t channels, size_t bins) 
{
    	fstream ff(flnm, std::fstream::in);
    	if(ff.fail()) 
	{	
		std::cerr << "\n ERROR: Failed to open INL file: " << flnm << endl; 
		return false;
	}

    	//Read the header from the file
    	//The format of the header seems to be [TDC-THESERIAL-inl_corr]
    	UInt_t TDCSerial;
    	ff.ignore(10, '-');
    	ff >> std::hex >> TDCSerial >> std::dec;
    	ff.ignore(1000, '\n');

	auto retval = m_mTDCs.insert(std::make_pair(TDCSerial, std::make_unique<MpdDeviceParameters>(channels, bins)));
	if(retval.second == false)
	{
            	cerr << "[MpdTofDigitProducer::Impl::setINLFromFile] TDC " << hex << TDCSerial<< dec << " already exist in the map."  << endl;
            	ff.close(); 
		return false;
	}

    	unsigned int lines_num = 0;
    	while(!ff.eof()) 
	{
		size_t chanId, binId;
            	string line; char dummy;

            	std::getline(ff, line, '\n');
            	if(ff.eof()) {return false;} 

		if(line == "") {continue;}
            	istringstream ss(line);

            	ss >> chanId >> dummy;
            	if(dummy != '=') 
		{
			cerr << "[MpdTofDigitProducer::Impl::setINLFromFile] Wrong INL file format." << endl; 
			ff.close(); 
			return false;
		}

		if(chanId > channels) 
		{
			cerr << "[MpdTofDigitProducer::Impl::setINLFromFile] Wrong channel in in the INL file." << endl; 
			ff.close(); 
			return false;
		}

		Double_t	value;
		while(ss.tellg() != -1) 
		{
                    	if(binId > bins) 
			{
                      		cerr << "[MpdTofDigitProducer::Impl::setINLFromFile] INL File contains too many bins in channel." << endl;
                            	ff.close(); 
				return false;
                    	}

			if(ss.peek()==',') {ss.ignore();}

                	ss >> value; 

			retval.first->second->SetINL(chanId, binId, value);

			binId++;		
            	}

            	if(binId != bins) 
		{
                    cout << "Warning: wrong number of bins in the INL file for channel " << chanId << " (" << binId << ")" << endl;
            	}

            	lines_num++;
    	}

    	if(lines_num != channels) 
	{
            cout << "Warning: wrong number of lines in the INL file (" << lines_num << endl;
	}

    	cout << "Tof: INL for TDC " << hex << TDCSerial << dec << " loaded succesfully from INL file." << endl;

return true;
}
//------------------------------------------------------------------------------------------------------------------------
void 		MpdTofDigitProducer::Impl::DumpMapping(const char* comment, std::ostream& os) const
{
	if(comment != nullptr) os<<comment;

	for(auto& entry : m_mTDCs)
	{
		os<<"\n serial="<<hex<<entry.first<<dec;
		entry.second->DumpMapping("", os);
	}
}
//------------------------------------------------------------------------------------------------------------------------
void 		MpdTofDigitProducer::Impl::DumpINLs(const char* comment, std::ostream& os) const
{
	if(comment != nullptr) os<<comment;

	for(auto& entry : m_mTDCs)
	{
		os<<"\n serial="<<hex<<entry.first<<dec;
		entry.second->DumpINLs("\n ", os);
	}
}
//------------------------------------------------------------------------------------------------------------------------

