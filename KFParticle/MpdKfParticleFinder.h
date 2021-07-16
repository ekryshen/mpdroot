#ifndef MPDKFPARTICLEFINDER_h
#define MPDKFPARTICLEFINDER_h

#define HomogeneousField
#define DO_TPCCATRACKER_EFF_PERFORMANCE

// ROOT includes 
#include <TClonesArray.h>
#include <TH1.h>
#include <TH2.h>

// FairRoot includes 
#include <FairTask.h>
#include <FairRunAna.h>
#include <FairMCEventHeader.h>
#include <FairFileSource.h>

// MPD includes 
#include <MpdTpcKalmanTrack.h>
#include <MpdKalmanFilter.h>
#include <MpdMCTrack.h>
#include <MpdVertex.h>
#include <MpdMiniDstFileSource.h>
#include <MpdMiniTrack.h>
#include <MpdMiniMcEvent.h>
#include <MpdMiniTrackCovMatrix.h>

// KFParticle includes
#include <KFPTrackVector.h>
#include <KFParticleTopoReconstructor.h>
#include <KFTopoPerformance.h>

// C++ includes 
#include <vector>
#include <map>

using namespace std;
using namespace TMath;

class ExtendedMpdTpcKalmanTrack : public MpdTpcKalmanTrack {
public:

    ExtendedMpdTpcKalmanTrack() : fIdx(-1) {
        ;
    }

    virtual ~ExtendedMpdTpcKalmanTrack() {
        ;
    }

    Int_t fIdx;

    ClassDef(ExtendedMpdTpcKalmanTrack, 0)
};

class TrackCuts : public TNamed {
public:
    // Setting cut values by default that correspond to the TPC kinem. acceptance ...

    TrackCuts() : fPtMin(.15), fPtMax(2.), fAbsEta(1.1) {
        ;
    }

    virtual ~TrackCuts() {
        ;
    }

    void SetPtMin(Double_t value) {
        fPtMin = value;
    }

    void SetPtMax(Double_t value) {
        fPtMax = value;
    }

    void SetAbsEta(Double_t value) {
        fAbsEta = value;
    }

    Double_t GetPtMin() {
        return fPtMin;
    }

    Double_t GetPtMax() {
        return fPtMax;
    }

    Double_t GetAbsEta() {
        return fAbsEta;
    }

private:
    Double_t fPtMin;
    Double_t fPtMax;

    Double_t fAbsEta;

    ClassDef(TrackCuts, 0)
};

class MpdKfParticleFinder : public FairTask {
public:

    MpdKfParticleFinder() {
        ;
    }
    MpdKfParticleFinder(TString);

    virtual ~MpdKfParticleFinder() {
        if (fTopoReconstructor)
            delete fTopoReconstructor;

        if (isUseCuts && fCuts)
            delete fCuts;
    }

    /// Init
    virtual InitStatus Init();

    /// Execute
    virtual void Exec(Option_t* option);

    /// Finish
    virtual void Finish();

    // User defined setters ...

    void SetPidHypo(Int_t pdg1, Int_t pdg2) {

        if (pdg1 < 0 || pdg2 > 0) {
            cout << "PDG1 is assumed to be positive, and PDG2 - negative!" << endl;
            throw;
        }

        fPidHypo.first = pdg1;
        fPidHypo.second = pdg2;
    }

    void SetUseCuts(Bool_t flag) {
        isUseCuts = flag;
    }

    void SetUseKFPerformance(Bool_t flag) {
        isPerformance = flag;
    }

    TrackCuts* GetTrackCuts() {
        if (!fCuts)
            return nullptr;
        else
            return fCuts;
    }

private:
    FairRootManager* ioman;

    pair <Int_t, Int_t> fPidHypo;
    vector <Int_t> fDecPdgs; // pdg codes of decays we are interested in ...

    TString fInput;
    vector <TString> fInputs;
    map <TString, UInt_t> fIds;

    TFile* fOutFile;

    TClonesArray* fTpcTracks;
    TClonesArray* fMcTracks;
    TClonesArray* fPrimaryVertices;
    FairMCEventHeader* fMCHeader;
    FairEventHeader* fHeader;

    // MiniDst
    TClonesArray* fMiniTracks;
    TClonesArray* fMiniCovMatrices;
    TClonesArray* fMCEvent;

    KFParticleTopoReconstructor* fTopoReconstructor;
    KFTopoPerformance* fPerformance;
    TDatabasePDG *fPdgDB;

    // Output arrays ...
    TClonesArray* fSecondaries;
    TClonesArray* fKFVertices;

    UInt_t evCounter = 0;

    // Cut flags ...
    Bool_t isUseCuts;
    TrackCuts* fCuts;

    Bool_t isMini;
    Bool_t isPerformance;

private:
    void ReadList(TString);

    Bool_t checkBranchStatus(TString dst) {
        // Checking branch status ...
        TChain* ch = new TChain(isMini ? "MiniDst" : "mpdsim");
        
        ch->Add(dst.Data());
        Bool_t status = ch->GetBranchStatus(isMini ? "Track" : "TpcKalmanTrack");

        delete ch;
        return status;
    }

    void Mini2Kalman(TClonesArray*);

    void ProcessDst();

    void FillKFPTrackVector(vector <ExtendedMpdTpcKalmanTrack>, KFPTrackVector*, vector <Int_t>);
    vector <KFMCTrack> FillKFMCTrack();
    Bool_t isInAcceptance(MpdMCTrack*);

    // Possible track params. to be used:
    void dcaToBeamline(MpdTpcKalmanTrack, MpdTpcKalmanTrack&); // orig. params written to output ...
    void tpcInnerShell(MpdTpcKalmanTrack, MpdTpcKalmanTrack&); // propagated track params. to  TPC inner shell ...
    void lastHit(MpdTpcKalmanTrack, MpdTpcKalmanTrack&); // track params. calculated at the last hit ...

    vector <Double_t> DoErrorPropagationToXYZPxPyPz(TMatrixDSym* cov, TMatrixD* params, MpdTpcKalmanTrack*);

    ClassDef(MpdKfParticleFinder, 0)
};

#endif