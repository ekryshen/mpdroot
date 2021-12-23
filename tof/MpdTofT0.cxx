//------------------------------------------------------------------------------------------------------------------------
/// \class MpdTofT0
/// 
/// \brief 
/// \author Viktor Baryshnikov (Faculty of Physics, MSU, Moscow)
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------
#include <assert.h>
#include <climits>


#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TEfficiency.h>
#include <TRandom2.h>
#include <TClonesArray.h>
#include <TStopwatch.h>

#include "MpdMCTrack.h"
#include "FairLogger.h"

#include "MpdTofUtils.h"
#include "MpdTofPoint.h"
#include "MpdTofHit.h"
#include "MpdTofGeoUtils.h"
#include "FairMCEventHeader.h"

#include "MpdTpcKalmanTrack.h"
//#include "MpdTofHit.h"
#include "MpdTofMatchingData.h"

#include "MpdTofT0.h"

using namespace std;

ClassImp(MpdTofT0)
ClassImp(MpdTofT0Data)
//------------------------------------------------------------------------------------------------------------------------
LState::LState(size_t& size)
{
	if(size > 15)
	{
		size = 15; // pow(3, 15) = 14,348,907 combination (~100 Mbyte)
		cout<<"\n -W: LState ctor: reduce size.";
	}

	fSize = size;
	fState.resize(size, 0); // init by 0
}
//------------------------------------------------------------------------------------------------------------------------
LState::~LState()
{
	
}
//------------------------------------------------------------------------------------------------------------------------
size_t		LState::Size() const
{
return pow(fBase, fSize);
}
//------------------------------------------------------------------------------------------------------------------------
void		LState::Next()
{
	fState[0]++;
	_checkAndShift(0);
}
//------------------------------------------------------------------------------------------------------------------------
void		LState::_checkAndShift(size_t index)
{
	if(fState[index] == fBase)
	{
		fState[index] = 0;
		fState[index+1]++;
	 	_checkAndShift(index+1);
	}
}
//------------------------------------------------------------------------------------------------------------------------
void 		LState::Print(const char* comment, ostream& os)const
{
	if(comment != nullptr) os<<comment;
	os<<" [";
	for(const auto& entry : fState)
	{
	 	os<<( (int)entry);
	}
	os<<"]";
}
//------------------------------------------------------------------------------------------------------------------------
MpdTofT0::MpdTofT0(const char *name, Bool_t useMCdata, Int_t verbose, const char *flnm)
  : FairTask(name, verbose), fUseMCData(useMCdata), fFlnm(flnm)
{
	fResult = new MpdTofT0Data("TOFTevent");
	fDoTest = (flnm != nullptr);
	pRandom = new TRandom2;

	if(fDoTest)
	{
		pTimer = new TStopwatch;

		Add(hNmatch = new TH1D("MpdTofT0_hNmatch", ";N tof matching per event;entries",  1000, -0.5, 999.5));
		Add(hTeveN = new TH2D("MpdTofT0_hTeveN", ";matching mupltiplicity; Tevent, ns",  1000, -0.5, 999.5, 1000, -5., 5.));
		Add(hTcalcN = new TH2D("MpdTofT0_hTcalcN", ";matching mupltiplicity; calc. time, s",  1000, -0.5, 999.5, 1000, 0., 1000.));
		Add(hTofmultB = new TH2D("MpdTofT0_hTofmultB", ";impact parameter, fm; matching mupltiplicity", 1000, 0., 20., 1000, -0.5, 999.5));

		Add(effT0_10 = new TEfficiency("MpdTofT0_effT0_10", "T0 < 10 ps;impact parameter, fm;efficiency",  20, 0., 20.));
		Add(effT0_15 = new TEfficiency("MpdTofT0_effT0_15", "T0 < 15 ps;impact parameter, fm;efficiency",  20, 0., 20.));
		Add(effT0_20 = new TEfficiency("MpdTofT0_effT0_20", "T0 < 20 ps;impact parameter, fm;efficiency",  20, 0., 20.));
		Add(effT0_50 = new TEfficiency("MpdTofT0_effT0_50", "T0 < 50 ps;impact parameter, fm;efficiency",  20, 0., 20.));
		Add(effT0_100 = new TEfficiency("MpdTofT0_effT0_100", "T0 < 100 ps;impact parameter, fm;efficiency",  20, 0., 20.));
	}
}
//------------------------------------------------------------------------------------------------------------------------
MpdTofT0::~MpdTofT0()
{
	delete pRandom;
	delete pTimer;
	fList.Delete();
}
//------------------------------------------------------------------------------------------------------------------------
InitStatus 		MpdTofT0::Init()
{
    		aTofMatchings = (TClonesArray*) FairRootManager::Instance()->GetObject("TOFMatching");
    		aTPCkfTracks = 	(TClonesArray*) FairRootManager::Instance()->GetObject("TpcKalmanTrack");
assert(aTofMatchings);
assert(aTPCkfTracks);

		FairRootManager::Instance()->Register("TOFTevent", "Tof", fResult, kFALSE); // if kTRUE, branch will be saved to the tree file

	if(fUseMCData)
	{
   		aMcTracks = 	(TClonesArray*) FairRootManager::Instance()->GetObject("MCTrack");  
		fEventHeader = (FairMCEventHeader*) FairRootManager::Instance()->GetObject("MCEventHeader."); 
assert(aMcTracks);
assert(fEventHeader);
	}

	LOG(INFO)<<"[MpdTofT0::Init] Initialization finished succesfully.";

return kSUCCESS;
}
//------------------------------------------------------------------------------------------------------------------------
void 			MpdTofT0::Exec(Option_t *option)
{
	// cleanup containers
	mmMatchData.clear();
	fResult->reset();

	struct {size_t Nmult = {};double impact = {}, Zvert = {};} 	MCdata;

	if(fUseMCData)
	{
		MCdata.Nmult = fEventHeader->GetNPrim();
		MCdata.impact = fEventHeader->GetB(); //[fm]
		MCdata.Zvert = fEventHeader->GetZ(); // [cm]

//		if(10. < MCdata.Zvert || MCdata.Zvert < -10.) return; 						// pass vertex Z inside[-10,10] cm <<<---------CUT!!!
	}

	if(pTimer) pTimer->Start();

	size_t TofMatchingNmb = aTofMatchings->GetEntriesFast();
	if(fDoTest) hNmatch->Fill(TofMatchingNmb);

	for(size_t i = 0; i < TofMatchingNmb; i++) // cycle by Tof matchings
	{
		auto match = (MpdTofMatchingData*) aTofMatchings->At(i);
		int KFtid = match->GetKFTrackIndex();
		auto pKfTrack = (MpdTpcKalmanTrack*) aTPCkfTracks->UncheckedAt(KFtid);	


		auto track = (MpdMCTrack*) aMcTracks->UncheckedAt(pKfTrack->GetTrackID());
		if(track->GetMotherId() != -1) continue; 								// pass primary tracks ONLY <<<---------CUT!!!

		mmMatchData.insert(make_pair(pKfTrack->Momentum3().Pt(), make_pair(match, pKfTrack)));
	}

	map<size_t, pair<size_t, TmatchMini[fRangeSize]> > mMatchMini; // < divisionID, <size, array> >

	// fill divisions by matchings(fRangeSize=15 matchings per divizion)
	size_t nMatchings = 0, divisionID = 0;
	for(auto it = mmMatchData.begin(), itEnd = mmMatchData.end(); it != itEnd; it++, nMatchings++) // cycle by of matchings, decrement by Pt
	{
		double P = it->second.second->Momentum();		
		double L = it->second.first->GetTrackLength();
		double Ttof = it->second.first->GetTime();

		if(nMatchings < fRangeSize) // add matching to current divizion
		{
			auto& data = mMatchMini[divisionID];
			data.first++;
			(data.second)[nMatchings] = TmatchMini(P, L, Ttof, it);
		}
		else // begin fill new divizion
		{
			divisionID++;
			nMatchings = 0;

			auto& data = mMatchMini[divisionID];
			data.first = 1;
			(data.second)[nMatchings] = TmatchMini(P, L, Ttof, it);
		}
	}

	double TevNom = 0., TevDenom = 0.;
	for(auto& divizion : mMatchMini) // cycle by divizions
	{	
		size_t size = divizion.second.first; // size <= fRangeSize

		if(size < 4) continue;
		auto result = Estimate(divizion.second.second, size);

		double weight = (size -2) / result.Chi2; //  1/(Chi2/NDF)

		TevNom += weight *  result.Tevent; // mean with weight
		TevDenom += weight;
	}

	if(pTimer)
	{
		pTimer->Stop();
		pTimer->Print("m");
	}

	if(TevDenom > 0.) // successfuly, update result
	{
		fResult->Tevent = TevNom / TevDenom;
		fResult->Chi2 = TevDenom; // sum of 1/(Chi2/NDF) weights
		fResult->nMatchings = mmMatchData.size();
	}

	if(fDoTest)
	{
		hTeveN->Fill(fResult->nMatchings, fResult->Tevent);
		hTcalcN->Fill(fResult->nMatchings, pTimer->CpuTime());

		if(fUseMCData)
		{
			double t0 = fabs(fResult->Tevent);
			effT0_10->Fill( fResult->Chi2 != 0. && t0 < 0.010, MCdata.impact);
			effT0_15->Fill( fResult->Chi2 != 0. && t0 < 0.015, MCdata.impact);
			effT0_20->Fill( fResult->Chi2 != 0. && t0 < 0.020, MCdata.impact);
			effT0_50->Fill( fResult->Chi2 != 0. && t0 < 0.050, MCdata.impact);
			effT0_100->Fill(fResult->Chi2 != 0. && t0 < 0.100, MCdata.impact);

//cout<<"\n t0= "<<fResult->Tevent<<", ch2="<<fResult->Chi2<<" b="<<MCdata.impact;

			hTofmultB->Fill(MCdata.impact, fResult->nMatchings);
		}		
	}

	LOG(INFO)<<"[MpdTofT0::Exec] matching nmb= "<<fResult->nMatchings<<", T event= "<<fResult->Tevent<<", weight= "<<fResult->Chi2;
}
//------------------------------------------------------------------------------------------------------------------------
MpdTofT0Data 		MpdTofT0::Estimate(T_matchMini* data, size_t size)
{
const static double c_inv = 1./(299792458. *1.e+2 *1.e-9); // [ns/cm]

	LState 		hypothesis(size); 
	MpdTofT0Data 	result;
	double 		minChi2 = 1.e+50; // big value

	for(size_t state = 0, stateSize = hypothesis.Size(); state < stateSize; state++) // cycle by hypothesis state
	{
		// calc. Tevent (numerator and denominator) for current hypothesis.
		double numerator = 0., denominator = 0.;

		for(size_t matching = 0; matching < size; matching++)
		{
			double P = data[matching].P;
			double L = data[matching].L;
			double Ttof = data[matching].Ttof;
			double m = hypothesis.Mass(matching);
			double Tsigma = hypothesis.TSigma(matching);

			double Texp = data[matching].Texp = sqrt(P*P + m*m) / P * c_inv * L; 
			double Weight = data[matching].Weight = 1./(fTofSigma*fTofSigma + Tsigma*Tsigma);

			numerator += Weight *(Ttof - Texp);
			denominator += Weight; 
		}

assert(denominator != 0.);	

		double Tevent = numerator / denominator;

		// calc. Chi2 for current hypothesis.
		double Chi2 = 0.;
		for(size_t matching = 0; matching < size; matching++)
		{
			double Weight = data[matching].Weight;
			double Ttof = data[matching].Ttof;
			double Texp = data[matching].Texp;

			// remove current matching deposit from Tevent
			double Tev = (numerator - Weight * (Ttof - Texp)) / (denominator - Weight);

			Chi2 += (Ttof - Tevent - Texp)*(Ttof - Tevent - Texp) * (Weight*Weight);
		}

		Chi2 /= size;

//cout<<"\n size="<<size<<" Chi2="<<Chi2;

		if(Chi2 < minChi2) // new best estimate of Tevent
		{
			result.Tevent = Tevent;
			result.Chi2 = Chi2;
			minChi2 = Chi2;

//cout<<"\n ----------->>>Local BEST minChi2="<<minChi2;
		}

		hypothesis.Next();

	} // cycle by hypothesis state

//cout<<"\n ------------------------------------------------------------------>>>> GLOBAL BEST minChi2="<<minChi2;

return result;
}
//------------------------------------------------------------------------------------------------------------------------
void 			MpdTofT0::Finish()
{
	if(fDoTest)
	{
		LOG(DEBUG2)<<"[MpdTofT0::Finish] Update  "<<fFlnm.Data()<<" file. ";
		auto ptr = gFile;
		TFile file(fFlnm.Data(), "RECREATE");
		fList.Write();
		file.Close();
		gFile = ptr;
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------
void 			MpdTofT0::SetSeed(UInt_t seed)
{
	pRandom->SetSeed(seed);
}
//--------------------------------------------------------------------------------------------------------------------------------------
template<typename T>
void			MpdTofT0::Add(T *hist)
{ 
assert(hist->InheritsFrom("TH1") || !std::string(hist->ClassName()).compare("TEfficiency"));

	hist->SetDirectory(nullptr); fList.Add(hist);
}
//------------------------------------------------------------------------------------------------------------------------


