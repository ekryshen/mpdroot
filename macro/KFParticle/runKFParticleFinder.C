#include <Rtypes.h>

R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "macro/mpd/mpdloadlibs.C"

void runKFParticleFinder(TString in = "", TString out = "", Int_t nStartEvent = 0, Int_t nEvents = 0) {
    FairRunAna* fRun = new FairRunAna();
    FairRootFileSink* sink = new FairRootFileSink(out);
    fRun->SetSink(sink);

    MpdKalmanFilter *kalman = MpdKalmanFilter::Instance("KF");
    fRun->AddTask(kalman);

    MpdKfParticleFinder* pF = new MpdKfParticleFinder(in);
    // pF->SetUseCuts(kTRUE);
    // pF->SetUseKFPerformance(kTRUE);
    
    // Setters to specify what a two-particle decay we are looking for ... 
    // pdg1, pdg2 to specify a decay. \Lambda^{0} is set by default
    // Available decays right now are: \Lambda^{0} -> p + \pion^{-} and \Kaon^{0}_{s} -> \pion^{+} + \pion^{-}
    // pF->SetPidHypo(+211, -211);
    
    fRun->AddTask(pF);
    fRun->Init();

    // Redefine established cuts if necessary ...
    //    TrackCuts* cutsGot = partFinder->GetTrackCuts();
    //    if (cutsGot) {
    //        cutsGot->SetPtMin(.5);
    //        cutsGot->SetPtMax(1.8);
    //        cutsGot->SetAbsEta(1.);
    //    }

    fRun->Run(nStartEvent, nStartEvent + nEvents);
    
    delete pF;
}