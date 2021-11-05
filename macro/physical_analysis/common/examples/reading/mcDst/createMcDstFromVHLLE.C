#include <Rtypes.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TClonesArray.h>
#include <TChain.h>

using namespace std;

R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "macro/mpd/mpdloadlibs.C"

/*
 A macro to be used for vHLLE input aiming at its converting to the McDst-format
 * @in@ means input vHLLE file
 * @out@ is considered as output file having *.mcDst.root extension
 */

void createMcDstFromVHLLE(TString in = "", TString out = "", TString branch = "treefin") {
    if (in.IsNull() || out.IsNull() || !out.Contains(".mcDst.")) {
        cout << "Specify correct input and output file names!" << endl;
        return;
    }
        
    // Prepare output tree with branches
    TFile* outFile = new TFile(out.Data(), "recreate");
    TTree* outTree = new TTree("McDst", "McDst");

    TClonesArray* events = new TClonesArray("McEvent");
    TClonesArray* particles = new TClonesArray("McParticle");

    outTree->Branch("Event", &events);
    outTree->Branch("Particle", &particles);

    // Open and read tree with vHLLE input
    TChain* fin = new TChain(branch.Data());
    fin->Add(in.Data());

    const Int_t dim = 15000;

    Float_t fX[dim];
    Float_t fY[dim];
    Float_t fZ[dim];
    Float_t fT[dim];

    Float_t fPx[dim];
    Float_t fPy[dim];
    Float_t fPz[dim];

    Float_t fE[dim];
    Int_t fID[dim];

    Int_t fMID[dim];
    Short_t fELE[dim];
    Short_t fBAR[dim];
    Short_t fSTR[dim];

    Int_t fNpart;

    fin->SetBranchAddress("x", fX);
    fin->SetBranchAddress("y", fY);
    fin->SetBranchAddress("z", fZ);
    fin->SetBranchAddress("t", fT);
    fin->SetBranchAddress("px", fPx);    
    fin->SetBranchAddress("py", fPy);
    fin->SetBranchAddress("pz", fPz);
    fin->SetBranchAddress("E", fE);
    fin->SetBranchAddress("id", fID);
    fin->SetBranchAddress("mid", fMID);
    fin->SetBranchAddress("ele", fELE);
    fin->SetBranchAddress("bar", fBAR);
    fin->SetBranchAddress("str", fSTR);
    fin->SetBranchAddress("npart", &fNpart);
   
    for (Int_t iEvent = 0; iEvent < fin->GetEntries(); iEvent++) {
        if (iEvent % 1000 == 0)
            cout << "Event# " << iEvent << endl;
        // Clear arrays to be used for writing
        events->Delete();
        particles->Delete();

        // Clear arrays to be used for reading
        for (Int_t iTrack = 0; iTrack < dim; iTrack++) {
            fX[iTrack] = fY[iTrack] = fZ[iTrack] = fT[iTrack] = 0.;
            fPx[iTrack] = fPy[iTrack] = fPz[iTrack] = fE[iTrack] = 0.;
            fID[iTrack] = fMID[iTrack] = fELE[iTrack] = fBAR[iTrack] = fSTR[iTrack] = fNpart = 0;
        }

        fin->GetEntry(iEvent);

        McEvent* event = new ((*events)[events->GetEntriesFast()]) McEvent();
        event->SetEventNr(iEvent);
        event->SetB(-1);

        for (Int_t iTrack = 0; iTrack < fNpart; iTrack++) {
            McParticle* particle = new ((*particles)[particles->GetEntriesFast()]) McParticle();

            particle->SetX(fX[iTrack]);
            particle->SetY(fY[iTrack]);
            particle->SetZ(fZ[iTrack]);

            particle->SetPx(fPx[iTrack]);
            particle->SetPy(fPy[iTrack]);
            particle->SetPz(fPz[iTrack]);

            particle->SetE(fE[iTrack]);
            particle->SetT(fT[iTrack]);

            particle->setPdg(fID[iTrack]);

        }
        outTree->Fill();
    }

    outTree->Write();
    outFile->Close();
}