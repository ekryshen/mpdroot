#include <TMath.h>
#include <TNamed.h>
#include <TFile.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <MpdTpcKalmanTrack.h>
#include <MpdTpcHit.h>
#include <MpdKalmanHit.h>
#include <TVector2.h>

#ifndef MPDFEMTOSHAREQUALITYPAIRCUT_H
#define MPDFEMTOSHAREQUALITYPAIRCUT_H 1

using namespace std;
using namespace TMath;

class MpdFemtoShareQualityPairCut : public TNamed {
public:

    MpdFemtoShareQualityPairCut() {;}
    MpdFemtoShareQualityPairCut(TFile*, Float_t, Float_t);
    MpdFemtoShareQualityPairCut(TChain*, TClonesArray*, TClonesArray*);

    virtual ~MpdFemtoShareQualityPairCut();

    // Getters

    UInt_t GetNPairsPassed() {
        return fNPairsPassed;
    }

    UInt_t GetNPairsFailed() {
        return fNPairsFailed;
    }
    
    Bool_t GetZeroSharing() {
        return fZeroSharing;
    }
    
    Float_t GetQualityMax() {
        return fShareQualityMax;
    }
    
    Float_t GetSharingMax() {
        return fShareFractionMax;
    }

    // Setters
    void SetShareQualityMax(Float_t val) {
        fShareQualityMax = val;
    }

    void SetShareFractionMax(Float_t val) {
        fShareFractionMax = val;
    }
    
    void SetZeroSharing(Bool_t val) {
        fZeroSharing = val;
    }
    
    void SetQualityMax(Float_t val) {
        fShareQualityMax = val;
    }

    void SetSharingMax(Float_t val) {
        fShareFractionMax = val;
    }
    
    void CheckTwoTrackEffects() {};
    Float_t Quality(Int_t, Int_t);
    Float_t Quality(Int_t, Int_t, UInt_t);
    Float_t Sharing(Int_t, Int_t);
    void MapOfSplittedTracks(TClonesArray*, map <Int_t, Int_t>&, map <Int_t, Int_t>::iterator&, Int_t);
    Bool_t Splitting(map <Int_t, Int_t>&, Int_t);
    void SplittingInefficiency(map <Int_t, Int_t>&) {}; 

private:

    inline Int_t GetHitIndex(TClonesArray* fHits, Int_t hit) {
        return ((MpdKalmanHit*) (fHits->UncheckedAt(hit)))->GetIndex(0);
    };

    inline Int_t GetHitLayer(TClonesArray* fHits, Int_t hit) {
        return ((MpdKalmanHit*) (fHits->UncheckedAt(hit)))->GetLayer();
    };

    TFile* fFile;
    TChain* fChain;

    TClonesArray* fTracks;
    TClonesArray* fMcTracks;
    TClonesArray* fHits_1st;
    TClonesArray* fHits_2nd;

    UInt_t fNPairsPassed; // Number of pairs consideered that passed the cut 
    UInt_t fNPairsFailed; // Number of pairs consideered that failed the cut

    Float_t fShareQualityMax; // Maximum allowed pair quality
    Float_t fShareFractionMax; // Maximum allowed share fraction
    
    Bool_t fZeroSharing;
     
    ClassDef(MpdFemtoShareQualityPairCut, 1)
};

#endif
