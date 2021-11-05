/* Macro reads DST file produced by macro reco.C */

#include <Rtypes.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TStopwatch.h>

R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "macro/mpd/mpdloadlibs.C"

void readDST(TString in = "") {
    if (in.IsNull()) {
        cout << "Please, provide an input DST-file!" << endl;
        return;
    }
          
    TStopwatch timer;
    timer.Start();

    TChain *dstTree = new TChain("mpdsim");
    dstTree->Add(in.Data());

    // Activate branches
    MpdEvent *event = nullptr;
    dstTree->SetBranchAddress("MPDEvent.", &event);
    TClonesArray *fMCTracks = nullptr;
    dstTree->SetBranchAddress("MCTrack", &fMCTracks);

    Int_t events = dstTree->GetEntries();
    cout << " Number of events in DST file = " << events << endl;

    for (Int_t i = 0; i < events; i++) {
        dstTree->GetEntry(i);

        Int_t Ntracks = event->GetGlobalTracks()->GetEntriesFast();

        cout << " Number of tracks = " << Ntracks << endl;
        for (Int_t iTrack = 0; iTrack < Ntracks; iTrack++) {
            MpdTrack* track = (MpdTrack*) event->GetGlobalTracks()->UncheckedAt(iTrack);
            
            //track->Dump();
            /*cout << "Track id = " << pDSTtrack->GetID()
             << "  Pt = "    << pDSTtrack->GetPt()
             << "  Theta = " << pDSTtrack->GetTheta()
             << "  Phi = "   << pDSTtrack->GetPhi()
             << "  Px = "  << pDSTtrack->GetPx()
             << "  Py = "  << pDSTtrack->GetPy()
             << "  Pz = "  << pDSTtrack->GetPz()
             << "  Eta = "  << pDSTtrack->GetEta()
             << endl;
             */
            /* See mpddata/MpdTrack.h for more methods */

        } // track loop
    } // event loop

    timer.Print();

    cout << " Test passed" << endl;
    cout << " All ok " << endl;
    exit(0);
}
