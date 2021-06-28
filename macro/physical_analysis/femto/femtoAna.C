#include <iostream>
#include <TChain.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <fstream>
#include <TH1.h>

using namespace std;

R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "macro/mpd/mpdloadlibs.C"

void femtoAna(TString inFileDST = "",
        TString outFile = "TEST_OLD.root", Int_t nStartEvent = 0, Int_t nEvents = 250) {

    Float_t Qinv = 0.20;
    MpdFemtoHistos* histos = new MpdFemtoHistos(Qinv, outFile);
    
    MpdFemto* femto = new MpdFemto(inFileDST, histos);

    femto->SetStartEvent(nStartEvent);
    femto->SetEvNumToRead(nEvents);
    femto->SetSourceSize(3.); // Source size in fm;
    femto->SetPdgCode(211); // pi+ are considered;
    femto->SetEtaCuts(-1.0, 1.0); // Define a range over pseudorapidity: val1 < eta < val2;
    femto->SetPtCuts(0.15, 1.5); // Define a range over Pt: val1 < Pt < val2;
    femto->SetNumMixedEvents(10); // Number of events to be mixed
    femto->SetKtCuts(0.15, 1.5); // Define a range over Kt: val1 < Kt < val2;
    femto->SetQinv(Qinv);
    femto->SetMagField(0.5); // Define magnitude of the mag.field
    femto->SetRadTpc(1.0); // Define average TPC rad. where split. && merg. effects are studied
    femto->SetMinNoHits(10);
    femto->SetZeroSharing(kTRUE); // If true, pairs with the number of shared hits equal to 0, are put in the histo _hSharing.  
    femto->SetSplittingCut(kTRUE);
    femto->SetSharingCut(kTRUE);
    femto->SetQualityMax(-0.4);
    femto->SetSharingMax(0.);

    // A method to perform 1D-analysis
    femto->MakeCFs_1D();

    // A method to perform DEtaDPhi-analysis
    // femto->DeltaEtaDeltaPhi();

    // A method to obtain Quality & Sharing cuts
    //femto->QualityAndSharing();

    delete femto;

    // Return parameters of the 1D-fit. [0] corresponds to R0, [1] -- to \lambda 
    Double_t* fitRes1D = histos->GetFitParams1D();
    cout << "R0 = " << fitRes1D[0] << " Lambda = " << fitRes1D[1] << endl;

    // Return parameters of the 3D-fit. [0] corresponds to Rout, [1] -- to Rside, [2] -- to Rlong
    //Double_t* fitRes3D = histos->GetFitParams3D();
    // cout << "Rout = " << fitRes3D[0] << " Rside = " << fitRes3D[1] << " Rlong = " << fitRes3D[2] << endl;  

    delete histos;
}
