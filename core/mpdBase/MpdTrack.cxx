// Author: Oleg Rogachevsky
// Update: 2009-09-17 18:43:28+0400
// Copyright: 2009 (C) MPD coll.
//
// Track container 

#ifndef ROOT_MpdTrack
#include "MpdTrack.h"
#include "FairRunAna.h"
#include "TMath.h"
#include "FairRunAna.h"
#include "FairField.h"
#include "Rtypes.h"
#endif



// -----   Default constructor ---------------------------------------
MpdTrack::MpdTrack():
  TObject(), fEdgeCut(0), fPidProbElectron(0.), fPidProbPion(0.),fPidProbKaon(0.), fPidProbProton(0.), 
  fPidTPCProbElectron(0.), fPidTPCProbPion(0.),fPidTPCProbKaon(0.), fPidTPCProbProton(0.), 
  fPidTOFProbElectron(0.), fPidTOFProbPion(0.),fPidTOFProbKaon(0.), fPidTOFProbProton(0.),      
  fTofBeta(0.), fTofMass2(0.), fdEdXTPC(0.), fTofFlag(0) ,fHitMap(0) ,
  fSharedHitMap(0),fNSigmaElectron(-100.), fNSigmaPion(-100.), fNSigmaKaon(-100.), fNSigmaProton(-100.), 
  fTofHitIndex(-1)  
{}
// -------------------------------------------------------------------
Float_t MpdTrack::GetPx() const
{
  return TMath::Abs(fPt)*TMath::Cos(fPhi); 
}
Float_t MpdTrack::GetPy() const
{
  return TMath::Abs(fPt)*TMath::Sin(fPhi);
}
Float_t MpdTrack::GetPz() const
{
  if ( TMath::Sin(fTheta) == 0.) return TMath::Sqrt(-1); // NaN
  return TMath::Abs(fPt)/TMath::Tan(fTheta);
}
Float_t MpdTrack::GetEta()const
{
  return -TMath::Log(TMath::Tan(0.5*fTheta));
}

// -------------------------------------------------------------------
ClassImp(MpdTrack);

MpdHelix MpdTrack::GetHelix() const {
	TVector3 mom(GetPx(),GetPy(),GetPz());
	TVector3 pos(GetFirstPointX(),GetFirstPointY(),GetFirstPointZ());
	Double_t charge = GetCharge();
	Double_t Bz = 0.5;
	if(FairRunAna::Instance()){
		if(FairRunAna::Instance()->GetField())
			Bz = FairRunAna::Instance()->GetField()->GetBz(0,0,0)*0.1;
	}
	return MpdHelix(mom,pos,charge,Bz);
}

Int_t MpdTrack::GetNSharedTpcHits() const {
	Int_t shared = 0;
	for(int i=0;i<53;i++)
		if(TESTBIT(fSharedHitMap,i))
			shared++;
	return shared;
}
