#ifndef MPDEMCHIT_H
#define MPDEMCHIT_H 1

#include "FairHit.h"

class MpdEmcHit : public FairHit {
public:

    /** Default constructor **/
    MpdEmcHit();

    /** Constructor with hit parameters (1)**/
    MpdEmcHit(Int_t detectorID, TVector3 pos, TVector3 dpos, Int_t refIndex, Int_t flag);

    /** Constructor with hit parameters (2) [not the flag]**/
    MpdEmcHit(Int_t detectorID, TVector3 pos, TVector3 dpos, Int_t refIndex);

    MpdEmcHit(UInt_t sec, UInt_t row, UInt_t supMod, UInt_t mod, Float_t e) {} //TEMPORARY FIX: was not implemented -> undefined reference
    MpdEmcHit(UInt_t sec, UInt_t row, UInt_t supMod, UInt_t mod, Float_t e, Float_t time);
    MpdEmcHit(UInt_t detID, Float_t e, Float_t time); 

    virtual ~MpdEmcHit();

    void Print(const Option_t* opt = 0) const;

    Int_t GetFlag() const {
        return fFlag;
    };

    Int_t GetSec() const {
        return fDetectorID/(1000*12);
    };

    Int_t GetMod() const {
       return fDetectorID - (fDetectorID/1000)*1000;
    };

    Int_t GetSupMod() const {
       return -1;	
    };

    Int_t GetRow() const {
        return fDetectorID/1000;;
    };

    Float_t GetE() const {
        return fE;
    };

    Float_t GetTime() const {
        return fTime;
    };

    Float_t GetRhoCenter() const {
        return sqrt(fX*fX+fY*fY);
    }

    Float_t GetZCenter() const {
        return fZ;
    };

    Float_t GetPhiCenter() const {
        return fPhiCenter;
    };

    Float_t GetThetaCenter() const {
        return fThetaCenter;
    };

    Int_t GetTrackId() const {
        return fTrackID;
    };

    Int_t GetPdg() const {
        return fPDG;
    };

    Int_t GetNumTracks() const {
        return fNumTracks;
    };

    Int_t GetDetectorID() const {
	return fDetectorID; 
    };

    Float_t GetX() const {
        return fX;
    };

    Float_t GetY() const {
        return fY;
    };

    Float_t GetZ() const {
        return fZ;
    };

    void SetFlag(Int_t flag) {
        fFlag = flag;
    };

    void SetEnergy(Float_t e) {
        fE = e;
    };

    void SetTime(Float_t time) {
	fTime = time;     
    };

    void IncreaseEnergy(Float_t e) {
        fE += e;
    };

    void IncreaseEnergyTime(Float_t timeEnergy) {
        fTime += timeEnergy;
    };

    void SetTrackId(Int_t id) {
        fTrackID = id;
    };

    void SetPdg(Int_t pdg) {
        fPDG = pdg;
    };

    void SetNumTracks(Int_t n) {
        fNumTracks = n;
    };

    void SetPhiCenter(Float_t phi) {
        fPhiCenter = phi;
    };

    void SetZCenter(Float_t z) {
        fZ = z;
    };

    void SetThetaCenter(Float_t theta) {
        fThetaCenter = theta;
    };

    void SetX(Float_t x) {
       fX = x; 
    };

    void SetY(Float_t y) {
       fY = y; 
    };

    void SetZ(Float_t z) {
       fZ = z; 
    };


protected:

    Float_t fE; //energy
    Float_t fTime; // hit mean time
    Int_t fDetectorID; // detector id of each hit 
    UInt_t fNumTracks; // number of tracks, which have contribution in module
    Float_t fPhiCenter; // phi-angle of the center of module
    Float_t fThetaCenter; // theta-angle of the center of module

    Int_t fTrackID; // -1 if more than one track in module
    Int_t fFlag; // Flag for general purposes [TDC, event tagging...]
    Int_t fPDG; // code of particle if only one track presented in module 

    ClassDef(MpdEmcHit, 1)

};


#endif
