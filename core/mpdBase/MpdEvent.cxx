// Author: Emelyanov D.
// Update: 2009-09-17 17:29:46+0400
// Copyright: 2009 (C) MPD coll.
//

#ifndef ROOT_MpdEvent
#include "MpdEvent.h"
#endif

// -----   Default constructor ---------------------------------------
MpdEvent::MpdEvent():
  TNamed("MpdEvent","Global"),
  fEventInfoNofPrimaryTracks(0),
  fEventInfoNofGlobalTracks(0),
  fPrimaryTracks(new TClonesArray("MpdTrack")),
  fGlobalTracks(new TClonesArray("MpdTrack"))
{}

// -------------------------------------------------------------------
MpdEvent::~MpdEvent()
{
  delete fPrimaryTracks;
  delete fGlobalTracks;
}

void MpdEvent::Reset(){
  fEventInfoNofPrimaryTracks = 0;
  fEventInfoNofGlobalTracks = 0;
  fPrimaryTracks->Delete();
  fGlobalTracks->Delete();
}

// -------------------------------------------------------------------
MpdTrack *MpdEvent::AddPrimaryTrack() 
{
  return new((*fPrimaryTracks)[fEventInfoNofPrimaryTracks++]) MpdTrack();
}

// -------------------------------------------------------------------
MpdTrack *MpdEvent::AddGlobalTrack()
{
  return new((*fGlobalTracks)[fEventInfoNofGlobalTracks++]) MpdTrack();
}
// -------------------------------------------------------------------

MpdEvent& MpdEvent::operator =(const MpdEvent& event) {
	if(&event==this){
		return *this;
	}
    fRunInfoRunId = event.fRunInfoRunId;
    fRunInfo_ProductionVersion = event.fRunInfo_ProductionVersion;
    fRunInfo_CenterOfMassEnergy = event.fRunInfo_CenterOfMassEnergy;
    fRunInfo_BeamMassNumber[0] = event.fRunInfo_BeamMassNumber[0];
    fRunInfo_BeamMassNumber[1] = event.fRunInfo_BeamMassNumber[1];
    fRunInfoMagneticFieldZ = event.fRunInfoMagneticFieldZ;
    fEventInfo_TriggerMask = event.fEventInfo_TriggerMask;
    fEventInfoNofGlobalTracks = event.fEventInfoNofGlobalTracks;
    fEventInfoNofPrimaryTracks = event.fEventInfoNofPrimaryTracks;
    fEventInfo_NofPositiveTracks = event.fEventInfo_NofPositiveTracks;
    fEventInfo_NofNegativeTracks = event.fEventInfo_NofNegativeTracks;

    // PrimaryVertices

    PrimaryVerticesX = event.PrimaryVerticesX;
    PrimaryVerticesY = event.PrimaryVerticesY;
    PrimaryVerticesZ = event.PrimaryVerticesZ;
    PrimaryVerticesXerr = event.PrimaryVerticesXerr;
    PrimaryVerticesYerr = event.PrimaryVerticesYerr;
    PrimaryVerticesZerr = event.PrimaryVerticesZerr;
    PrimaryVertices_SumTrackPt = event.PrimaryVertices_SumTrackPt;
    PrimaryVertices_MeanDip = event.PrimaryVertices_MeanDip;
    PrimaryVerticesChi2 = event.PrimaryVerticesChi2;
    PrimaryVertices_RefMultNeg = event.PrimaryVertices_RefMultNeg;
    PrimaryVertices_RefMultPos = event.PrimaryVertices_RefMultPos;
    fPrimaryTracks->Clear();
    fGlobalTracks->Clear();
    for(int i=0;i<event.fPrimaryTracks->GetEntriesFast();i++){
    	MpdTrack *from = (MpdTrack*)event.fPrimaryTracks->UncheckedAt(i);
    	MpdTrack *to = (MpdTrack*)fPrimaryTracks->ConstructedAt(i);
    	*to = *from;
    }
    for(int i=0;i<event.fGlobalTracks->GetEntriesFast();i++){
    	MpdTrack *from = (MpdTrack*)event.fGlobalTracks->UncheckedAt(i);
    	MpdTrack *to = (MpdTrack*)fGlobalTracks->ConstructedAt(i);
    	*to = *from;
    }
    return *this;
}

// -------------------------------------------------------------------

ClassImp(MpdEvent);


