#include "MpdFemto.h"
#include "MpdFemtoContainer.h"

//--------------------------------------------------------------------------

MpdFemto::MpdFemto() {
    fPDG = 0;
    fMass = 0;
    fFilename = "";
    fEtaCutLow = 0.;
    fEtaCutUp = 0.;
    fPtCutLow = 0.;
    fPtCutUp = 0.;
    fKtCutLow = 0.;
    fKtCutUp = 0.;
    fSourceSize = 0.;
    fMixedEvents = 0;
}

//--------------------------------------------------------------------------

MpdFemto::MpdFemto(const Char_t* fname) :
fFilename(""),
fPartTable(NULL),
fMassPart(NULL),
fDstTree(NULL),
fMpdEvent(NULL),
fMcTracks(NULL),
fRecoTracks(NULL),
fFemtoContainerReco(NULL),
fFemtoContainerMc(NULL),
fMpdTrackReco(NULL),
fMpdTrackMc(NULL) {
    TString path = getenv("VMCWORKDIR");
    TString path2pdg = path + "/input/pdg_table.txt";

    fPartTable = new TDatabasePDG();
    fPartTable->ReadPDGTable(path2pdg.Data());

    fFilename = fname;
    fDstTree = new TChain("cbmsim");
    fDstTree->Add(fFilename);

    Int_t events = fDstTree->GetEntries();
    cout << "Number of events in " << fname << " is " << events << endl;
    // cout << endl;

    fDstTree->SetBranchAddress("MPDEvent.", &fMpdEvent);
    fDstTree->SetBranchAddress("MCTrack", &fMcTracks);
    fFemtoContainerReco = new TClonesArray("MpdFemtoContainer");
    fFemtoContainerMc = new TClonesArray("MpdFemtoContainer");

    fHisto = new MpdFemtoHistos();
}

//--------------------------------------------------------------------------

MpdFemto::MpdFemto(const Char_t* fname, MpdFemtoHistos* _h) :
fFilename(""),
fPartTable(NULL),
fMassPart(NULL),
fDstTree(NULL),
fMpdEvent(NULL),
fMcTracks(NULL),
fRecoTracks(NULL),
fFemtoContainerReco(NULL),
fFemtoContainerMc(NULL),
fMpdTrackReco(NULL),
fMpdTrackMc(NULL) {
    TString path = getenv("VMCWORKDIR");
    TString path2pdg = path + "/input/pdg_table.txt";

    fPartTable = new TDatabasePDG();
    fPartTable->ReadPDGTable(path2pdg.Data());

    fFilename = fname;
    fDstTree = new TChain("cbmsim");
    fDstTree->Add(fFilename);

    Int_t events = fDstTree->GetEntries();
    cout << "Number of events in " << fname << " is " << events << endl;
    // cout << endl;

    fDstTree->SetBranchAddress("MPDEvent.", &fMpdEvent);
    fDstTree->SetBranchAddress("MCTrack", &fMcTracks);
    fFemtoContainerReco = new TClonesArray("MpdFemtoContainer");
    fFemtoContainerMc = new TClonesArray("MpdFemtoContainer");

    fHisto = _h;
}
//--------------------------------------------------------------------------

MpdFemto::~MpdFemto() {
    delete fDstTree;
    delete fFemtoContainerReco;
    delete fFemtoContainerMc;
    delete fPartTable;
    // delete fHisto;
}

//--------------------------------------------------------------------------

void MpdFemto::ReadEvent(Int_t evNum) {
    
    fDstTree->GetEntry(evNum);

    fRecoTracks = fMpdEvent->GetGlobalTracks();
    // cout << "Number of RECO-tracks: " << fRecoTracks->GetEntriesFast() << endl;
    // cout << "Number of MC-tracks: " << fMcTracks->GetEntriesFast() << endl;

    for (Int_t iRecoTrack = 0; iRecoTrack < fRecoTracks->GetEntriesFast(); iRecoTrack++) {

        fMpdTrackReco = (MpdTrack*) fRecoTracks->UncheckedAt(iRecoTrack);
        Int_t trackID = fMpdTrackReco->GetID();

        fMpdTrackMc = (FairMCTrack*) fMcTracks->UncheckedAt(trackID);

        if (fMpdTrackMc->GetPdgCode() != fPDG)
            continue;

        if (fMpdTrackReco->GetEta() < fEtaCutLow || fMpdTrackReco->GetEta() > fEtaCutUp)
            continue;

        if (Abs(fMpdTrackReco->GetPt()) < fPtCutLow || Abs(fMpdTrackReco->GetPt()) > fPtCutUp)
            continue;

        Float_t Px_reco = fMpdTrackReco->GetPx();
        Float_t Py_reco = fMpdTrackReco->GetPy();
        Float_t Pz_reco = fMpdTrackReco->GetPz();
        Float_t E_reco = Sqrt(Px_reco * Px_reco + Py_reco * Py_reco + Pz_reco * Pz_reco + fMass * fMass);

        Float_t Px_sim = fMpdTrackMc->GetPx();
        Float_t Py_sim = fMpdTrackMc->GetPy();
        Float_t Pz_sim = fMpdTrackMc->GetPz();
        Float_t E_sim = Sqrt(Px_sim * Px_sim + Py_sim * Py_sim + Pz_sim * Pz_sim + fMass * fMass);

        Float_t x = gRandom->Gaus(0, fSourceSize * Sqrt(2.0));
        Float_t y = x;
        Float_t z = y;

        TLorentzVector MOM_RECO(Px_reco, Py_reco, Pz_reco, E_reco);
        TLorentzVector MOM_MC(Px_sim, Py_sim, Pz_sim, E_sim);
        TLorentzVector COORD(x, y, z, 0.);

        new((*fFemtoContainerReco)[fFemtoContainerReco->GetEntriesFast()])
                MpdFemtoContainer(evNum, MOM_RECO, COORD);

        //        fFemtoContainerReco->AddLast(new MpdFemtoContainer(evNum, MOM_RECO, COORD));
        //        fFemtoContainerMc->AddLast(new MpdFemtoContainer(evNum, MOM_MC, COORD));

        new((*fFemtoContainerMc)[fFemtoContainerMc->GetEntriesFast()])
                MpdFemtoContainer(evNum, MOM_MC, COORD);
    }
}

Float_t MpdFemto::EposFemtoWeightQS(TLorentzVector a, TLorentzVector b, TLorentzVector x, TLorentzVector y) {
    return 1. + Cos(5.068423 * ((a.X() - b.X()) * (x.X() - y.X())
            + (a.Y() - b.Y()) * (x.Y() - y.Y())
            + (a.Z() - b.Z()) * (x.Z() - y.Z())
            - (a.T() - b.T()) * (x.T() - y.T())
            ));
}

void MpdFemto::MakeCFs_1D() {

    cout << "1D-analysis started ... " << endl;
    cout << "Number of events to be mixed: " << fMixedEvents << endl;

    fMassPart = fPartTable->GetParticle(fPDG);
    fMass = fMassPart->Mass();

    for (Int_t iEvent = 0; iEvent < fEvNum; iEvent += fMixedEvents) {
        fFemtoContainerMc->Delete();
        fFemtoContainerReco->Delete();

        for (Int_t jEvent = iEvent; jEvent < iEvent + fMixedEvents; jEvent++)
            ReadEvent(jEvent);

        cout << "Event: " << iEvent << " of " << fDstTree->GetEntries() << " mix " << fFemtoContainerReco->GetEntriesFast() << " particles" << endl;

        for (Int_t iPart = 0; iPart < fFemtoContainerReco->GetEntriesFast(); iPart++)
            for (Int_t jPart = iPart + 1; jPart < fFemtoContainerReco->GetEntriesFast(); jPart++) {

                TLorentzVector mom_iPart_reco = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(iPart))->Get4Momentum();
                TLorentzVector mom_jPart_reco = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(jPart))->Get4Momentum();
                TLorentzVector mom_iPart_sim =  ((MpdFemtoContainer*) fFemtoContainerMc->UncheckedAt(iPart))->Get4Momentum();
                TLorentzVector mom_jPart_sim =  ((MpdFemtoContainer*) fFemtoContainerMc->UncheckedAt(jPart))->Get4Momentum();

                TLorentzVector coord_iPart_reco = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(iPart))->Get4Coordinate();
                TLorentzVector coord_jPart_reco = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(jPart))->Get4Coordinate();

                Double_t Qinv = EposFemtoQinv4vec(mom_iPart_reco, mom_jPart_reco);

                TVector3 vlong;
                vlong.SetXYZ(0., 0., (mom_iPart_reco.Z() + mom_jPart_reco.Z()) / (mom_iPart_reco.E() + mom_jPart_reco.E()));

                mom_iPart_reco.Boost(-vlong);
                mom_jPart_reco.Boost(-vlong);

                Double_t kt = 0.5 * Hypot(mom_iPart_reco.X() + mom_jPart_reco.X(), mom_iPart_reco.Y() + mom_jPart_reco.Y());
                if (kt < fKtCutLow || kt > fKtCutUp)
                    continue;

                Int_t iPart_evNum = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(iPart))->GetEventNumber();
                Int_t jPart_evNum = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(jPart))->GetEventNumber();

                if (iPart_evNum == jPart_evNum) {
                    Float_t wfemto = EposFemtoWeightQS(mom_iPart_sim, mom_jPart_sim, coord_iPart_reco, coord_jPart_reco);
                    fHisto->GetNominator()->Fill(Qinv, wfemto);
                    fHisto->GetNominatorBase()->Fill(Qinv);

                } else
                    fHisto->GetDenominator()->Fill(Qinv);
            }
    }
}

void MpdFemto::MakeCFs_3D() {
 
    cout << "3D-analysis started ... " << endl;
    cout << "Number of events to be mixed: " << fMixedEvents << endl;

    fMassPart = fPartTable->GetParticle(fPDG);
    fMass = fMassPart->Mass();

    for (Int_t iEvent = 0; iEvent < fEvNum; iEvent += fMixedEvents) {
        fFemtoContainerMc->Delete();
        fFemtoContainerReco->Delete();

        for (Int_t jEvent = iEvent; jEvent < iEvent + fMixedEvents; jEvent++)
            ReadEvent(jEvent);

        cout << "Event: " << iEvent << " of " << fDstTree->GetEntries() << " mix " << fFemtoContainerReco->GetEntriesFast() << " particles" << endl;

        for (Int_t iPart = 0; iPart < fFemtoContainerReco->GetEntriesFast(); iPart++)
            for (Int_t jPart = iPart + 1; jPart < fFemtoContainerReco->GetEntriesFast(); jPart++) {

                TLorentzVector mom_iPart_reco = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(iPart))->Get4Momentum();
                TLorentzVector mom_jPart_reco = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(jPart))->Get4Momentum();
                TLorentzVector mom_iPart_sim = ((MpdFemtoContainer*) fFemtoContainerMc->UncheckedAt(iPart))->Get4Momentum();
                TLorentzVector mom_jPart_sim = ((MpdFemtoContainer*) fFemtoContainerMc->UncheckedAt(jPart))->Get4Momentum();

                TLorentzVector coord_iPart_reco = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(iPart))->Get4Coordinate();
                TLorentzVector coord_jPart_reco = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(jPart))->Get4Coordinate();

                TVector3 vlong;
                vlong.SetXYZ(0., 0., (mom_iPart_reco.Z() + mom_jPart_reco.Z()) / (mom_iPart_reco.E() + mom_jPart_reco.E()));

                mom_iPart_reco.Boost(-vlong);
                mom_jPart_reco.Boost(-vlong);

                Double_t kt = 0.5 * Hypot(mom_iPart_reco.X() + mom_jPart_reco.X(), mom_iPart_reco.Y() + mom_jPart_reco.Y());
                if (kt < fKtCutLow || kt > fKtCutUp)
                    continue;

                Int_t iPart_evNum = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(iPart))->GetEventNumber();
                Int_t jPart_evNum = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(jPart))->GetEventNumber();

                Float_t pXsum = mom_iPart_reco.X() + mom_jPart_reco.X();
                Float_t pYsum = mom_iPart_reco.Y() + mom_jPart_reco.Y();
                Float_t pXdif = mom_iPart_reco.X() - mom_jPart_reco.X();
                Float_t pYdif = mom_iPart_reco.Y() - mom_jPart_reco.Y();

                Float_t qOut = 0.5 * (pXsum * pXdif + pYsum * pYdif) / kt;
                Float_t qSide = (mom_iPart_reco.X() * mom_jPart_reco.Y() - mom_iPart_reco.Y() * mom_jPart_reco.X()) / kt;
                Float_t qLong = mom_iPart_reco.Z() - mom_jPart_reco.Z();
                
                if (qOut < 0. && qOut > fQinv && qSide < 0. && qSide > fQinv && qLong < 0. && qLong > fQinv)
                    continue;

                if (iPart_evNum == jPart_evNum) {
                    Float_t wfemto = EposFemtoWeightQS(mom_iPart_sim, mom_jPart_sim, coord_iPart_reco, coord_jPart_reco);
                    fHisto->GetNominator3D()->Fill(qOut, qSide, qLong, wfemto);
                } else
                    fHisto->GetDenominator3D()->Fill(qOut, qSide, qLong);
            }
    }
}