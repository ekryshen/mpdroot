#include <iostream>
#include <TNamed.h>
#include <TFile.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <MpdEvent.h>
#include <MpdMCTrack.h>
#include <FairMCPoint.h>
#include <MpdTrack.h>
#include <MpdTpcDigit.h>
#include <TMath.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <MpdFemtoHistos.h>
#include "MpdFemtoShareQualityPairCut.h"
#include <TVector3.h>

#ifndef MPDFEMTO_H
#define MPDFEMTO_H 1

using namespace std;
using namespace TMath;

class MpdFemto : public TNamed {
public:
    MpdFemto();
    MpdFemto(const Char_t*);
    MpdFemto(const Char_t*, MpdFemtoHistos*);
    MpdFemto(const Char_t*, const Char_t*, MpdFemtoHistos*);
    virtual ~MpdFemto();

    // Getters

    Int_t GetPdgCode() {
        return fPDG;
    }

    Int_t GetEntriesNum() {
        return fDstTree->GetEntries();
    }

    MpdFemtoHistos* GetHistos() {
        return fHisto;
    }

    Float_t GetQinv() {
        return fQinv;
    }
    
    Float_t GetMagField() {
        return fMagField;
    }
    
    Float_t GetRadTpc() {
        return fRadTPC;
    }
    
    Int_t GetStartEvent() {
        return fStartEvent;
    }

    // Setters

    void SetPdgCode(Int_t val) {
        fPDG = val;
    }

    void SetEtaCuts(Float_t low, Float_t up) {
        fEtaCutLow = low;
        fEtaCutUp = up;
    }

    void SetPtCuts(Float_t low, Float_t up) {
        fPtCutLow = low;
        fPtCutUp = up;
    }

    void SetKtCuts(Float_t low, Float_t up) {
        fKtCutLow = low;
        fKtCutUp = up;
    }

    void SetSourceSize(Float_t size) {
        fSourceSize = size;
    }

    void SetNumMixedEvents(Int_t num) {
        fMixedEvents = num;
    }

    void SetQinv(Float_t qinv) {
        fQinv = qinv;
    }

    void SetNbins(Int_t val) {
        fBins = val;
    }

    void SetUpLimit(Float_t xUp) {
        fxUp = xUp;
    }

    void SetEvNumToRead(Int_t val) {
        fEvNum = val;
    }
    
    void SetMagField(Float_t val) {
        fMagField = val;
    }
    
    void SetRadTpc(Float_t val) {
        fRadTPC = val;
    } 
    
    void SetStartEvent(Int_t val) {
        fStartEvent = val;
    }
    
    void SetSharingCut(Bool_t flag) {
        fSharingCut = flag;
    }
    
    void SetSplittingCut(Bool_t flag) {
        fSplittingCut = flag;
    }

    void SetMinNoHits(Int_t val) {
        fMinNoHits = val;
    }
    
    void SetZeroSharing(Bool_t val) {
        fCuts->SetZeroSharing(val);
    }
    
    void SetQualityMax(Float_t val) {
        fCuts->SetQualityMax(val);
    }

    void SetSharingMax(Float_t val) {
        fCuts->SetSharingMax(val);
    }
    
    void SetDeltaEtaDeltaPhi(Bool_t val) {
        fDeltaEtaDeltaPhi = val;
    }
       
    void MakeCFs_1D();
    void MakeCFs_3D() {};
    void DeltaEtaDeltaPhi();
    
    void QualityAndSharing();
    void EffSplitting() {;}
    
    inline void DebugInfo() {
        cout << " Service information: " << endl;
        cout << " fPDG = " << fPDG << endl;
        cout << " fMass = " << fMass << endl;
        cout << " fQinv = " << fQinv << endl;
        cout << " fKtCutLow = " << fKtCutLow << endl;
        cout << " fKtCutUp = " << fKtCutUp << endl;
        cout << " fEtaCutLow = " << fEtaCutLow << endl;
        cout << " fEtaCutUp = " << fEtaCutUp << endl;
        cout << " fPtCutLow = " << fPtCutLow << endl;
        cout << " fPtCutUp = " << fPtCutUp << endl;
        cout << " fSourceSize = " << fSourceSize << endl;
        cout << " fMixedEvents = " << fMixedEvents << endl;
        cout << " fMcTracks, size = " << fMcTracks->GetEntriesFast() << endl;
        // cout << " fRecoTracks, size = " << fRecoTracks->GetEntriesFast() << endl;
        // cout << " fFemtoContainerReco, size = " << fFemtoContainerReco->GetEntriesFast() << endl;
        // cout << " fFemtoContainerMc, size = " << fFemtoContainerMc->GetEntriesFast() << endl;
    }


private:
    MpdFemtoHistos* fHisto;
    MpdFemtoShareQualityPairCut* fCuts;
    
    Int_t fMinNoHits;
    
    Bool_t fSharingCut;
    Bool_t fSplittingCut;
    Bool_t fDeltaEtaDeltaPhi;
    // Bool_t fZeroSharing;
    
    void ReadEvent(Int_t, const Char_t*);
    Bool_t Preselection(MpdTpcKalmanTrack*, MpdMCTrack*);
       
    Int_t fPDG;
    Float_t fMass;
    Float_t fCharge;
    Float_t fMagField;
    Float_t fRadTPC;
    const Char_t* fFilenameEve;
    const Char_t* fFilenameDST;
    TDatabasePDG* fPartTable;
    TParticlePDG* fParticle;

    Float_t fQinv;
    Float_t fKtCutLow;
    Float_t fKtCutUp;
    Float_t fEtaCutLow;
    Float_t fEtaCutUp;
    Float_t fPtCutLow;
    Float_t fPtCutUp;
    Float_t fSourceSize;
    Int_t fMixedEvents;

    Int_t fBins;
    Float_t fxUp;

    Int_t fEvNum;
    Int_t fStartEvent;

    TChain* fDstTree;
    TChain* fEveTree;
    MpdEvent* fMpdEvent;

    TClonesArray* fMcTracks;
    TClonesArray* fMcPoints;
    TClonesArray* fRecoTracks;
    TClonesArray* fFemtoContainerReco;
    TClonesArray* fFemtoContainerMc;

    MpdTrack* fMpdTrackReco;
    MpdMCTrack* fMpdTrackMc;
        
    TClonesArray* fTracksTPC;
    TClonesArray* fDigitsTPC;
    // MpdTpcKalmanTrack* fKalmanTrack;
    // TClonesArray* fHitsTPC;
    
   // MpdTpcKalmanTrack* fTrack;
   // MpdKalmanHit* fKalmanHit;

    inline Float_t EposFemtoQinv4vec(TLorentzVector first, TLorentzVector second) {
        return Abs((first - second).M());
    }

    Float_t EposFemtoWeightQS(TLorentzVector, TLorentzVector, TLorentzVector, TLorentzVector);

    ClassDef(MpdFemto, 1)
};

#endif
