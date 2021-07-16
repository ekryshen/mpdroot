#include <Rtypes.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>

R__ADD_INCLUDE_PATH($VMCWORKDIR)

void readKFParticleOutput(TString in = "", Int_t nEvents = 0) {
    TChain* ch = new TChain("mpdsim");
    ch->Add(in.Data());

    TClonesArray* primVertices = nullptr;
    TClonesArray* recoSecondaries = nullptr;
    FairEventHeader* eventHeader = nullptr;

    ch->SetBranchAddress("EventHeader.", &eventHeader);
    ch->SetBranchAddress("Vertex", &primVertices);
    ch->SetBranchAddress("Secondaries", &recoSecondaries);

    // Example of possible output histograms ...

    enum Proj {
        X, Y, Z
    };

    TH1F* hVpX = new TH1F(Form("%d", Proj::X), "V_{p}(X); V_{p}(X) [cm]; ", 100, -2., +2.);
    TH1F* hVpY = new TH1F(Form("%d", Proj::Y), "V_{p}(Y); V_{p}(Y) [cm]; ", 100, -2., +2.);
    TH1F* hVpZ = new TH1F(Form("%d", Proj::Z), "V_{p}(Z); V_{p}(Z) [cm]; ", 50, -70., +70.);

    TH1F* hMassWithoutCuts = new TH1F("hMassWithoutCuts", "hMassWithoutCuts", 75, 1.07, 1.22);
    TH1F* hMassWithCuts = new TH1F("hMassWithCuts", "hMassWithCuts", 75, 1.07, 1.22);

    TH2F* hPtEta = new TH2F("Pt vs. #eta", "Pt vs. #eta; #eta; Pt [GeV/c]", 100, -1.2, +1.2, 100, 0.2, 2.5);

    for (Int_t iEvent = 0; iEvent < (nEvents ? nEvents : ch->GetEntries()); iEvent++) {
        ch->GetEntry(iEvent);

        if (iEvent % 1000 == 0)
            cout << "Processing event# " << iEvent << endl;

        // Getting information on reco file and event used ...
        Int_t fileId = eventHeader->GetInputFileId(); // id of file in the list processed by runKFParticle.C
        Int_t evNum = eventHeader->GetMCEntryNumber(); // current event that corresponds to the same event in reco output 

        // Getting reconstructed primary vertex in event ...
        // In general case, here should be a loop since we are able to reconstruct 
        // more than one primary vertex ...
        // One vertes is set by default in macro/KFParticle/runKFParticleFinder.C 
        // It means that primVertices->GetEntriesFast() is equal to 1
        for (Int_t iVertex = 0; iVertex < primVertices->GetEntriesFast(); iVertex++) {
            KFVertex* Vp = (KFVertex*) primVertices->UncheckedAt(iVertex);

            // To get a list of all possible accessors to written data, please, look at KFParticle/KFParticle/KFVertex.h

            Float_t* covMatrix = Vp->CovarianceMatrix();

            // Getting components of reconstructed vertex ...
            Double_t X = Vp->GetX();
            Double_t Y = Vp->GetY();
            Double_t Z = Vp->GetZ();

            // X = Y = Z = 0. means that primary vertex has not been reconstructed in event ...
            if (!X || !Y || !Z)
                continue;

            hVpX->Fill(X);
            hVpY->Fill(Y);
            hVpZ->Fill(Z);

            // Int_t nContributors = Vp->GetNContributors();
            // cout << nContributors << endl;
        }

        // Getting reconstructed secondaries ...
        // Pdg hypo of a particle we are working with and other cuts to be used ...
        const Int_t pdg = 3122;
        const Double_t PtMin = 0.5;
        const Double_t PtMax = 2.;
        const Double_t AbsEta = 1.1;

        for (Int_t iSecondary = 0; iSecondary < recoSecondaries->GetEntriesFast(); iSecondary++) {
            KFParticle* particle = (KFParticle*) recoSecondaries->UncheckedAt(iSecondary);

            // To get a list of all possible accessors to written data, please, look at KFParticle/KFParticle/KFParticle.h and KFParticle/KFParticle/KFParticleBase.h

            Int_t pdgHypo = particle->GetPDG();
            if (pdgHypo != pdg)
                continue;

            Double_t mass = particle->GetMass();

            hMassWithoutCuts->Fill(mass);

            Double_t Px = particle->GetPx();
            Double_t Py = particle->GetPy();
            Double_t Pz = particle->GetPz();

            Double_t pT = particle->GetPt();
            Double_t eta = particle->GetEta();

            if (pT < PtMin || pT > PtMax || TMath::Abs(eta) > AbsEta)
                continue;

            hMassWithCuts->Fill(mass);
            hPtEta->Fill(eta, pT);
            
            // It is possible to get indices of reconstructed tracks in standard / MiniDst. 
            // Doing it for standard dst, one can get source Monte Carlo track to check 
            // whether its mother is a desirable decaying particle we are looking for ...
            // (daughter id given by particle finder --> idx of reconstructed track --> 
            // --> idx of corresponding MC track --> idx of mother particle)
            // To do it, one has to link standard dst to the macro in an usual way 
            // fileId and evNum are available to get corresponding file, event and arrays with reconstructed tracks 
                       
            const Double_t massThresh = 0.025;
            if (TMath::Abs(mass) < massThresh) {
                Int_t nDaugh = particle->NDaughters();
                vector <Int_t> daughters = particle->DaughterIds();

                for (Int_t iDaugh = 0; iDaugh < nDaugh; iDaugh++)
                    cout << daughters[iDaugh] << " ";
                cout << endl;
            }
        }
    }

    TCanvas* c = new TCanvas("c", "c", 1200, 800);
    c->Divide(3, 2);

    c->cd(1);
    hVpX->Draw();

    c->cd(2);
    hVpY->Draw();

    c->cd(3);
    hVpZ->Draw();

    c->cd(4);
    hMassWithoutCuts->Draw("PE1X0");

    c->cd(5);
    hMassWithCuts->Draw("PE1X0");

    c->cd(6);
    hPtEta->Draw("colz");

    c->SaveAs("test.pdf");

    delete ch;
}