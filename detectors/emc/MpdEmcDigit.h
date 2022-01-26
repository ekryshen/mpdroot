//--------------------------------------------------------------------
//
// Description:
//      MPD EMC Digit
//
//
// Author List:
//      Alexander Zinchenko LHEP, JINR, Dubna - 8-May-2016
//
//--------------------------------------------------------------------

#ifndef MPDEMCDIGIT_H
#define MPDEMCDIGIT_H 1

#include "FairHit.h"
#include <map>

class MpdEmcDigit : public FairHit 
{
 public:
  enum Shifts {kSideS = 0, kChannelS = 1, kRowS = 7, kSectorS = 13};
  enum Masks {kSideM = 1, kChannelM = 63, kRowM = 63, kSectorM = 15};

 public:

    /** Default constructor **/
    MpdEmcDigit();

    /** Constructor with hit parameters (1)**/
    MpdEmcDigit(Int_t detectorID, TVector3 pos, TVector3 dpos, Int_t refIndex, Int_t flag);

    /** Constructor with hit parameters (2) [not the flag]**/
    MpdEmcDigit(Int_t detectorID, TVector3 pos, TVector3 dpos, Int_t refIndex);

    MpdEmcDigit(UInt_t side, UInt_t sector, UInt_t row, UInt_t channel);

    /** Destructor **/
    virtual ~MpdEmcDigit();

    void Print(const Option_t* opt = 0) const;
    void IncreaseEnergy(Float_t e, Int_t trId);
    UInt_t SetDetId(UInt_t side, UInt_t sector, UInt_t module, UInt_t channel);
    UInt_t Side() const { return fDetectorID % 2; }
    UInt_t Channel() const { return (fDetectorID >> kChannelS) & kChannelM; }
    UInt_t Row() const { return (fDetectorID >> kRowS) & kRowM; }
    UInt_t Sector() const { return (fDetectorID >> kSectorS) & kSectorM; }

    /** Accessors **/
    Int_t GetFlag() const { return fFlag; }
    Int_t GetChanZId() const { return fChanZId; }
    Int_t GetChanPhiId() const { return fChanPhiId; }
    Float_t GetE() const { return fE; }
    Float_t GetZcenter() const { return fZCenter; }
    Float_t GetPhiCenter() const { return fPhiCenter; }
    Int_t GetTrackId() const { return fTrackID; }
    Int_t GetPdg() const { return fPDG; }
    Int_t GetNumTracks() const { return fNumTracks; }
    //std::map<Int_t,Float_t>& GetContrib() { return fContrib; }
    std::map<Int_t,Float_t> GetContrib() { return fContrib; }

    /** Modifiers **/
    void SetFlag(Int_t flag) { fFlag = flag; }
    void SetEnergy(Float_t e) { fE = e; }
    void IncreaseEnergy(Float_t e) { fE += e; }
    void SetTrackId(Int_t id) { fTrackID = id; }
    void SetPdg(Int_t pdg) { fPDG = pdg; }
    void SetNumTracks(Int_t n) { fNumTracks = n; }
    void SetZCenter(Float_t z) { fZCenter = z; }
    void SetPhiCenter(Float_t phi) { fPhiCenter = phi; }
    void SetChanZId(UInt_t zId) { fChanZId = zId; }
    void SetChanPhiId(UInt_t phiId) { fChanPhiId = phiId; }

protected:

    // base digi-part
    Int_t fChanZId; // channel number along Z
    Int_t fChanPhiId; // channel number along phi
    Float_t fE; //energy

    //extra reco-part
    Int_t fTrackID; // -1 if more than one track in module
    Int_t fFlag; // Flag for general purposes [TDC, event tagging...]
    Int_t fPDG; // code of particle if only one track presented in module 
    UInt_t fNumTracks; // number of tracks, which have contribution in module
    Float_t fZCenter; // z-coordinate of the center of module
    Float_t fPhiCenter; // phi-angle of the center of module
    std::map<Int_t,Float_t> fContrib; // energy deposit for each track ID

    ClassDef(MpdEmcDigit, 1);

};


#endif
