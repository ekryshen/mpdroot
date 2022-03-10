// Author: Emelyanov D.
// Update: Oleg Rogachevsky 2009-09-17 17:21:45+0400
// Copyright: 2009 (C) MPD coll.
//

#ifndef ROOT_MpdEvent
#define ROOT_MpdEvent
#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include "TNamed.h"
#include <TClonesArray.h>
#include "MpdTrack.h"

class MpdEvent : public TNamed {

 private:

    Int_t       fRunInfoRunId;
    Int_t       fRunInfo_ProductionVersion;
    Float_t     fRunInfo_CenterOfMassEnergy;
    Int_t       fRunInfo_BeamMassNumber[2];
    Double_t    fRunInfoMagneticFieldZ;
    UInt_t      fEventInfo_TriggerMask;
    Int_t       fEventInfoNofGlobalTracks;
    Int_t       fEventInfoNofPrimaryTracks;
    Int_t       fEventInfo_NofPositiveTracks;
    Int_t       fEventInfo_NofNegativeTracks;

    // PrimaryVertices

    Float_t     PrimaryVerticesX;
    Float_t     PrimaryVerticesY;
    Float_t     PrimaryVerticesZ;
    Float_t     PrimaryVerticesXerr;
    Float_t     PrimaryVerticesYerr;
    Float_t     PrimaryVerticesZerr;
    Float_t     PrimaryVertices_SumTrackPt;
    Float_t     PrimaryVertices_MeanDip;
    Float_t     PrimaryVerticesChi2;
    UShort_t    PrimaryVertices_RefMultNeg;
    UShort_t    PrimaryVertices_RefMultPos;

    TClonesArray *fPrimaryTracks; //-> array of only primary tracks
    TClonesArray *fGlobalTracks; //-> array of all tracks

 public:
 
  MpdEvent(); // Default constructor
  MpdEvent& operator=(const MpdEvent &event);
  virtual ~MpdEvent(); // Destructor

  void Reset();

    // SETTERS

    void SetRunInfoRunId(Int_t n) { fRunInfoRunId = n;}
    void SetRunInfo_ProductionVersion(Int_t n) { fRunInfo_ProductionVersion = n;}
    void SetRunInfo_CenterOfMassEnergy(Float_t n) {fRunInfo_CenterOfMassEnergy = n;}
    void SetRunInfo_BeamMassNumber(Int_t n, Int_t i) {fRunInfo_BeamMassNumber[i] = n;}
    void SetRunInfoMagneticFieldZ(Double_t n) {fRunInfoMagneticFieldZ = n;}
    void SetEventInfo_TriggerMask(UInt_t n) {fEventInfo_TriggerMask = n;}
    void SetEventInfoNofGlobalTracks(Int_t n) {fEventInfoNofGlobalTracks = n;}
    void SetEventInfoNofPrimaryTracks(Int_t n) {fEventInfoNofPrimaryTracks = n;}
    void SetEventInfo_NofPositiveTracks(Int_t n) {fEventInfo_NofPositiveTracks = n;}
    void SetEventInfo_NofNegativeTracks(Int_t n) {fEventInfo_NofNegativeTracks = n;}

    // PrimaryVertices

    void SetPrimaryVerticesX(Float_t n) {PrimaryVerticesX = n;}
    void SetPrimaryVerticesY(Float_t n) {PrimaryVerticesY = n;}
    void SetPrimaryVerticesZ(Float_t n) {PrimaryVerticesZ = n;}
    void SetPrimaryVerticesXerr(Float_t n) {PrimaryVerticesXerr = n;}
    void SetPrimaryVerticesYerr(Float_t n) {PrimaryVerticesYerr = n;}
    void SetPrimaryVerticesZerr(Float_t n) {PrimaryVerticesZerr = n;}
    void SetPrimaryVertices_SumTrackPt(Float_t n) {PrimaryVertices_SumTrackPt = n;}
    void SetPrimaryVertices_MeanDip(Float_t n) {PrimaryVertices_MeanDip = n;}
    void SetPrimaryVerticesChi2(Float_t n) {PrimaryVerticesChi2 = n;}
    void SetPrimaryVertices_RefMultNeg(UShort_t n) {PrimaryVertices_RefMultNeg = n;}
    void SetPrimaryVertices_RefMultPos(UShort_t n) {PrimaryVertices_RefMultPos = n;}

    // GETTERS

    Int_t   GetRunInfoRunId() {return fRunInfoRunId;}
    Int_t   GetRunInfo_ProductionVersion() {return fRunInfo_ProductionVersion;}
    Float_t     GetRunInfo_CenterOfMassEnergy() {return fRunInfo_CenterOfMassEnergy;}
    Int_t   GetRunInfo_BeamMassNumber(int n) {return fRunInfo_BeamMassNumber[n];}	// !!! n 0 or 1
    Double_t    GetRunInfoMagneticFieldZ() {return fRunInfoMagneticFieldZ;}
    UInt_t  GetEventInfo_TriggerMask() {return fEventInfo_TriggerMask;}
    Int_t   GetEventInfoNofGlobalTracks() {return fEventInfoNofGlobalTracks;}
    Int_t   GetEventInfoNofPrimaryTracks() {return fEventInfoNofPrimaryTracks;}
    Int_t   GetEventInfo_NofPositiveTracks() {return fEventInfo_NofPositiveTracks;}
    Int_t   GetEventInfo_NofNegativeTracks() {return fEventInfo_NofNegativeTracks;}
    
    TClonesArray *GetPrimaryTracks() {return fPrimaryTracks;} 
    TClonesArray *GetGlobalTracks() {return fGlobalTracks;} 


    // PrimaryVertices

    Float_t     GetPrimaryVerticesX() {return PrimaryVerticesX;}
    Float_t     GetPrimaryVerticesY() {return PrimaryVerticesY;}
    Float_t     GetPrimaryVerticesZ() {return PrimaryVerticesZ;}
    Float_t     GetPrimaryVerticesXerr() {return PrimaryVerticesXerr;}
    Float_t     GetPrimaryVerticesYerr() {return PrimaryVerticesYerr;}
    Float_t     GetPrimaryVerticesZerr() {return PrimaryVerticesZerr;}
    Float_t     GetPrimaryVertices_SumTrackPt() {return PrimaryVertices_SumTrackPt;}
    Float_t     GetPrimaryVertices_MeanDip() {return PrimaryVertices_MeanDip;}
    Float_t     GetPrimaryVerticesChi2() {return PrimaryVerticesChi2;}
    UShort_t    GetPrimaryVertices_RefMultNeg() {return PrimaryVertices_RefMultNeg;}
    UShort_t    GetPrimaryVertices_RefMultPos() {return PrimaryVertices_RefMultPos;}

    MpdTrack *AddPrimaryTrack();
    MpdTrack *AddGlobalTrack();

  ClassDef(MpdEvent,1) //MPDEVENT
};

#endif
