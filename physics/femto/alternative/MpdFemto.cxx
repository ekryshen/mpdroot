#include <Rtypes.h>

#include "MpdFemto.h"
#include "MpdFemtoContainer.h"
#include "MpdFemtoShareQualityPairCut.h"

//--------------------------------------------------------------------------

MpdFemto::MpdFemto() {
    fPDG = 0;
    fMass = 0;
    fFilenameDST = "";
    fFilenameEve = "";
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
fFilenameEve(""),
fFilenameDST(""),
fPartTable(NULL),
fParticle(NULL),
fDstTree(NULL),
fMpdEvent(NULL),
fMcTracks(NULL),
fRecoTracks(NULL),
fFemtoContainerReco(NULL),
fFemtoContainerMc(NULL),
fMpdTrackReco(NULL),
fMpdTrackMc(NULL),
fMass(0.),
fSplittingCut(kFALSE),
fSharingCut(kFALSE), 
fDeltaEtaDeltaPhi(kFALSE) {
    TString path = getenv("VMCWORKDIR");
    TString path2pdg = path + "/input/pdg_table.txt";

    fPartTable = new TDatabasePDG();
    fPartTable->ReadPDGTable(path2pdg.Data());

    fFilenameDST = fname;
    fDstTree = new TChain("mpdsim");
    fDstTree->Add(fFilenameDST);

    Int_t events = fDstTree->GetEntries();
    cout << "Number of events in " << fname << " is " << events << endl;
    // cout << endl;

    fDstTree->SetBranchAddress("MPDEvent.", &fMpdEvent);
    fDstTree->SetBranchAddress("MCTrack", &fMcTracks);
    // fDstTree->SetBranchAddress("TpcKalmanTrack", &fTracksTPC);

    fFemtoContainerReco = new TClonesArray("MpdFemtoContainer");
    fFemtoContainerMc = new TClonesArray("MpdFemtoContainer");

    fHisto = new MpdFemtoHistos();
    //fCuts = new MpdFemtoShareQualityPairCut(fDstTree, fTracksTPC);   
}

//--------------------------------------------------------------------------

MpdFemto::MpdFemto(const Char_t* fname, MpdFemtoHistos* _h) :
fFilenameEve(""),
fFilenameDST(""),
fPartTable(NULL),
fParticle(NULL),
fDstTree(NULL),
fMpdEvent(NULL),
fMcTracks(NULL),
fRecoTracks(NULL),
fFemtoContainerReco(NULL),
fFemtoContainerMc(NULL),
fMpdTrackReco(NULL),
fMpdTrackMc(NULL),
fTracksTPC(NULL),
fMass(0.),
fMinNoHits(0),
fSplittingCut(kFALSE),
fSharingCut(kFALSE), 
fDeltaEtaDeltaPhi(kFALSE) {
    TString path = getenv("VMCWORKDIR");
    TString path2pdg = path + "/input/pdg_table.txt";

    fPartTable = new TDatabasePDG();
    fPartTable->ReadPDGTable(path2pdg.Data());

    fFilenameDST = fname;
    fDstTree = new TChain("mpdsim");
    fDstTree->Add(fFilenameDST);

    Int_t events = fDstTree->GetEntries();
    cout << "Number of events in " << fname << " is " << events << endl;
    // cout << endl;

    fDstTree->SetBranchAddress("MPDEvent.", &fMpdEvent);
    fDstTree->SetBranchAddress("MCTrack", &fMcTracks);
    fDstTree->SetBranchAddress("TpcKalmanTrack", &fTracksTPC);

    fFemtoContainerReco = new TClonesArray("MpdFemtoContainer");
    fFemtoContainerMc = new TClonesArray("MpdFemtoContainer");

    fHisto = _h;
    fCuts = new MpdFemtoShareQualityPairCut(fDstTree, fTracksTPC, fMcTracks);
}

MpdFemto::MpdFemto(const Char_t* fnameDST, const Char_t* fnameEve, MpdFemtoHistos* _h) :
fFilenameEve(""),
fFilenameDST(""),
fPartTable(NULL),
fParticle(NULL),
fDstTree(NULL),
fMpdEvent(NULL),
fMcTracks(NULL),
fRecoTracks(NULL),
fFemtoContainerReco(NULL),
fFemtoContainerMc(NULL),
fMpdTrackReco(NULL),
fMpdTrackMc(NULL),
fTracksTPC(NULL),
fMass(0.),
fMinNoHits(0),
fMcPoints(NULL),
fDigitsTPC(NULL),
fSplittingCut(kFALSE),
fSharingCut(kFALSE),
fDeltaEtaDeltaPhi(kFALSE) {
    TString path = getenv("VMCWORKDIR");
    TString path2pdg = path + "/input/pdg_table.txt";

    fPartTable = new TDatabasePDG();
    fPartTable->ReadPDGTable(path2pdg.Data());

    fFilenameDST = fnameDST;
    fDstTree = new TChain("mpdsim");
    fDstTree->Add(fFilenameDST);

    fFilenameEve = fnameEve;
    fEveTree = new TChain("mpdsim");
    fEveTree->Add(fFilenameEve);

    Int_t eventsDST = fDstTree->GetEntries();
    cout << "Number of events in " << fnameDST << " is " << eventsDST << endl;
    // cout << endl;

    Int_t eventsEve = fEveTree->GetEntries();
    cout << "Number of events in " << fnameEve << " is " << eventsEve << endl;
    // cout << endl;

    fDstTree->SetBranchAddress("MPDEvent.", &fMpdEvent);
    fDstTree->SetBranchAddress("MCTrack", &fMcTracks);
    fDstTree->SetBranchAddress("TpcKalmanTrack", &fTracksTPC);

    fEveTree->SetBranchAddress("TpcPoint", &fMcPoints);
    fEveTree->SetBranchAddress("MpdTpcDigit", &fDigitsTPC);

    fFemtoContainerReco = new TClonesArray("MpdFemtoContainer");
    fFemtoContainerMc = new TClonesArray("MpdFemtoContainer");

    fHisto = _h;
    fCuts = new MpdFemtoShareQualityPairCut(fDstTree, fTracksTPC, fMcTracks);
}

//--------------------------------------------------------------------------

MpdFemto::~MpdFemto() {
    delete fDstTree;
    delete fFemtoContainerReco;
    delete fFemtoContainerMc;
    delete fPartTable;
    delete fCuts;
    // delete fHisto;
}

//--------------------------------------------------------------------------

void MpdFemto::ReadEvent(Int_t evNum, const Char_t* track) {
    fDstTree->GetEntry(evNum);

    map <Int_t, Int_t> nPoinsInTrack;
    map <Int_t, Int_t>::iterator it;

    if (track == "globTracks") {
        fRecoTracks = fMpdEvent->GetGlobalTracks();

        for (Int_t iRecoTrack = 0; iRecoTrack < fRecoTracks->GetEntriesFast(); iRecoTrack++) {

            fMpdTrackReco = (MpdTrack*) fRecoTracks->UncheckedAt(iRecoTrack);
            Int_t trackID = fMpdTrackReco->GetID();

            fMpdTrackMc = (MpdMCTrack*) fMcTracks->UncheckedAt(trackID);

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

            Float_t Phi = fMpdTrackReco->GetPhi();
            Float_t Theta = fMpdTrackReco->GetTheta();

            new((*fFemtoContainerReco)[fFemtoContainerReco->GetEntriesFast()])
                    MpdFemtoContainer(evNum, MOM_RECO, COORD, Phi, Theta, iRecoTrack);

            new((*fFemtoContainerMc)[fFemtoContainerMc->GetEntriesFast()])
                    MpdFemtoContainer(evNum, MOM_MC, COORD, 0.0, 0.0, trackID);
        }
    } 
    
    else if (track == "kalmanTracks") {

        if (fSplittingCut)
	  fCuts->MapOfSplittedTracks(fTracksTPC, nPoinsInTrack, it, fMinNoHits);

        for (Int_t iKalmanTrack = 0; iKalmanTrack < fTracksTPC->GetEntriesFast(); iKalmanTrack++) {
            MpdTpcKalmanTrack* tr1 = (MpdTpcKalmanTrack*) fTracksTPC->UncheckedAt(iKalmanTrack);
            Int_t trId1 = tr1->GetTrackID();
            Int_t nHits1 = tr1->GetNofHits();
	    //cout << fMinNoHits << endl;
	    // Nantes, 07July 2016
	    if (nHits1 < fMinNoHits)
	      continue;
	    
            if (fSplittingCut)
                if (fCuts->Splitting(nPoinsInTrack, trId1))
                    continue;

            MpdMCTrack* mc1 = (MpdMCTrack*) fMcTracks->UncheckedAt(trId1);

            if (!Preselection(tr1, mc1))
                continue;

            if (fSharingCut) {
                Bool_t isShared = kFALSE;
                for (Int_t jKalmanTrack = 0; jKalmanTrack < fTracksTPC->GetEntriesFast(); jKalmanTrack++) {
                    if (iKalmanTrack == jKalmanTrack)
                        continue;

                    MpdTpcKalmanTrack* tr2 = (MpdTpcKalmanTrack*) fTracksTPC->UncheckedAt(jKalmanTrack);
                    Int_t trId2 = tr2->GetTrackID();
                    Int_t nHits2 = tr2->GetNofHits();
//
                    if (nHits1 != nHits2)
                        continue;

                    MpdMCTrack* mc2 = (MpdMCTrack*) fMcTracks->UncheckedAt(trId2);

                    if (!Preselection(tr2, mc2))
                        continue;

                    Float_t s = fCuts->Sharing(iKalmanTrack, jKalmanTrack);
                    Float_t q = fCuts->Quality(iKalmanTrack, jKalmanTrack);
                    fHisto->GetQuality()->Fill(q);
                    fHisto->GetSharing()->Fill(s);
                    fHisto->GetQualityVsSharing()->Fill(q, s);
                    fHisto->GetQualityVsNhits()->Fill(q, nHits1);
                    if (q > -0.4) {
                        isShared = kTRUE;
                        break;
                    }
                }
                if (isShared)
                    continue;
            }

            Float_t Px_reco = tr1->Momentum3().X();
            Float_t Py_reco = tr1->Momentum3().Y();
            Float_t Pz_reco = tr1->Momentum3().Z();
            Float_t E_reco = Sqrt(tr1->Momentum() * tr1->Momentum() + fMass * fMass);

            Float_t Px_sim = mc1->GetPx();
            Float_t Py_sim = mc1->GetPy();
            Float_t Pz_sim = mc1->GetPz();
            Float_t E_sim = Sqrt(Px_sim * Px_sim + Py_sim * Py_sim + Pz_sim * Pz_sim + fMass * fMass);

            Float_t x = gRandom->Gaus(0, fSourceSize * Sqrt(2.0));
            Float_t y = x;
            Float_t z = y;

            TLorentzVector MOM_RECO(Px_reco, Py_reco, Pz_reco, E_reco);
            TLorentzVector MOM_MC(Px_sim, Py_sim, Pz_sim, E_sim);
            TLorentzVector COORD(x, y, z, 0.);

            Float_t Phi = tr1->Phi();
            Float_t Theta = tr1->Theta();

            new((*fFemtoContainerReco)[fFemtoContainerReco->GetEntriesFast()])
                    MpdFemtoContainer(evNum, MOM_RECO, COORD, Phi, Theta, trId1);

            new((*fFemtoContainerMc)[fFemtoContainerMc->GetEntriesFast()])
                    MpdFemtoContainer(evNum, MOM_MC, COORD, 0.0, 0.0, trId1);
        }

    } else {
        cout << "Define a type of tracks to be used! " << endl;
        return;
    }
    if (fSplittingCut) {
        for (map <Int_t, Int_t>::const_iterator it1 = nPoinsInTrack.begin(); it1 != nPoinsInTrack.end(); ++it1) {
            MpdMCTrack* track1 = (MpdMCTrack*) fMcTracks->UncheckedAt(it1->first);
            fHisto->GetNsplits()->Fill(it1->second);
            if (it1->second > 1)
                fHisto->GetPtSplit()->Fill(track1->GetPt());
            else
                fHisto->GetPtNoSplit()->Fill(track1->GetPt());
        }
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

    fParticle = fPartTable->GetParticle(fPDG);
    fMass = fParticle->Mass();

    for (Int_t iEvent = fStartEvent; iEvent < fEvNum + fStartEvent; iEvent += fMixedEvents) {
        fFemtoContainerMc->Delete();
        fFemtoContainerReco->Delete();

        if (iEvent > fDstTree->GetEntries() - 1)
            return;

        for (Int_t jEvent = iEvent; jEvent < iEvent + fMixedEvents; jEvent++)
            ReadEvent(jEvent, "kalmanTracks");

        cout << "Event: " << iEvent << " of " << fDstTree->GetEntries() << " mix " << fFemtoContainerReco->GetEntriesFast() << " particles" << endl;

       Float_t phid, etad;//, eta1, eta2, phi1, phi2, arg1, arg2;
        for (Int_t iPart = 0; iPart < fFemtoContainerReco->GetEntriesFast(); iPart++) {
            

//            if (fDeltaEtaDeltaPhi) {
//               
//                fHisto->GetEtaPhiStar()->Fill(eta1, phi1);
//            }
                       
            for (Int_t jPart = iPart + 1; jPart < fFemtoContainerReco->GetEntriesFast(); jPart++) {
                TLorentzVector mom_iPart_reco = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(iPart))->Get4Momentum();
                TLorentzVector mom_iPart_sim = ((MpdFemtoContainer*) fFemtoContainerMc->UncheckedAt(iPart))->Get4Momentum();
                TLorentzVector coord_iPart_reco = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(iPart))->Get4Coordinate();
                TLorentzVector mom_jPart_reco = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(jPart))->Get4Momentum();
                TLorentzVector mom_jPart_sim = ((MpdFemtoContainer*) fFemtoContainerMc->UncheckedAt(jPart))->Get4Momentum();
                TLorentzVector coord_jPart_reco = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(jPart))->Get4Coordinate();
                 
                if (fDeltaEtaDeltaPhi) {
                    Float_t pt1 = Sqrt(mom_iPart_reco.X() * mom_iPart_reco.X() + mom_iPart_reco.Y() * mom_iPart_reco.Y());
                    Float_t phi1 = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(iPart))->GetPhi();
                    Float_t theta1 = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(iPart))->GetTheta();
                    Float_t eta1 = -Log(Tan(theta1 / 2.));
                    Float_t arg1 = -0.3 * fMagField * fCharge * fRadTPC / (2. * pt1);
                    Float_t pt2 = Sqrt(mom_jPart_reco.X() * mom_jPart_reco.X() + mom_jPart_reco.Y() * mom_jPart_reco.Y());
                    Float_t phi2 = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(jPart))->GetPhi();
                    Float_t theta2 = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(jPart))->GetTheta();
                    Float_t eta2 = -Log(Tan(theta2 / 2.));
                    Float_t arg2 = -0.3 * fMagField * fCharge * fRadTPC / (2. * pt2);

                    phid = phi2 - phi1 + ASin(arg2) - ASin(arg1); // cout << phid << endl;
                    etad = eta2 - eta1;

                    while (phid > Pi())
                        phid -= 2. * Pi();
                    while (phid < -Pi())
                        phid += 2. * Pi();
                    
                   // cout << phid << endl;
                }

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
                    if (fDeltaEtaDeltaPhi) 
                        if (Qinv < 0.1)
                            fHisto->GetDeltaEtaDeltaPhiNomin()->Fill(phid, etad);
                } 
                
                else {
                    fHisto->GetDenominator()->Fill(Qinv);
                    if (fDeltaEtaDeltaPhi)
                        if (Qinv < 0.1)
                            fHisto->GetDeltaEtaDeltaPhiDenom()->Fill(phid, etad);
                }
            }
        }
    }
}

void MpdFemto::DeltaEtaDeltaPhi() {
    cout << "DeltaEtaDeltaPhi-analysis started ... " << endl;

    fParticle = fPartTable->GetParticle(fPDG);

    fMass = fParticle->Mass();
    fCharge = fParticle->Charge() / 3.; // Charge() returns the charge value in |e| / 3

    Int_t nPassedPairs = 0;
    Int_t nFailedPairs = 0;

    for (Int_t iEvent = fStartEvent; iEvent < fStartEvent + fEvNum; iEvent++) {
        fFemtoContainerMc->Delete();
        fFemtoContainerReco->Delete();

        if (iEvent > fDstTree->GetEntries() - 1)
            return;

        cout << "Event: " << iEvent << " of " << fDstTree->GetEntries() << " in the file" << endl;
        ReadEvent(iEvent, "kalmanTracks");

        for (Int_t iPart = 0; iPart < fFemtoContainerReco->GetEntriesFast(); iPart++) {
            TLorentzVector mom_iPart_reco = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(iPart))->Get4Momentum();
            //           Int_t trId1 = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(iPart))->GetTrackID();

            Float_t pt1 = Sqrt(mom_iPart_reco.X() * mom_iPart_reco.X() + mom_iPart_reco.Y() * mom_iPart_reco.Y());
            Float_t phi1 = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(iPart))->GetPhi();
            Float_t theta1 = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(iPart))->GetTheta();
            Float_t eta1 = -Log(Tan(theta1 / 2.));
            Float_t arg1 = -0.3 * fMagField * fCharge * fRadTPC / (2. * pt1);

            for (Int_t jPart = iPart + 1; jPart < fFemtoContainerReco->GetEntriesFast(); jPart++) {
                TLorentzVector mom_jPart_reco = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(jPart))->Get4Momentum();

                Float_t pt2 = Sqrt(mom_jPart_reco.X() * mom_jPart_reco.X() + mom_jPart_reco.Y() * mom_jPart_reco.Y());
                Float_t phi2 = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(jPart))->GetPhi();
                Float_t theta2 = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(jPart))->GetTheta();
                Float_t eta2 = -Log(Tan(theta2 / 2.));
                Float_t arg2 = -0.3 * fMagField * fCharge * fRadTPC / (2. * pt2);

                Float_t phid = phi2 - phi1 + ASin(arg2) - ASin(arg1);
                Float_t etad = eta2 - eta1;

                while (phid > Pi())
                    phid -= 2. * Pi();
                while (phid < -Pi())
                    phid += 2. * Pi();

                Double_t Qinv = EposFemtoQinv4vec(mom_iPart_reco, mom_jPart_reco);

                if (Qinv < fQinv)
                    fHisto->GetDeltaEtaDeltaPhi()->Fill(phid, etad);
            }
        }
        // cout << "nGoodPairs = " << nPassedPairs << " nFailedPairs = " << nFailedPairs<< endl;
        nPassedPairs = 0;
        nFailedPairs = 0;
    }
    fHisto->GetDeltaEtaDeltaPhi()->Scale(1. / fDstTree->GetEntries());
}

void MpdFemto::QualityAndSharing() {
    cout << "Define Quality and Sharing ... " << endl;

    fParticle = fPartTable->GetParticle(fPDG);

    fMass = fParticle->Mass();
    fCharge = fParticle->Charge() / 3.; // Charge() returns the charge value in |e| / 3

    for (Int_t iEvent = fStartEvent; iEvent < fStartEvent + fEvNum; iEvent++) {
        fFemtoContainerMc->Delete();
        fFemtoContainerReco->Delete();

        if (iEvent > fDstTree->GetEntries() - 1)
            return;

        cout << "Event: " << iEvent << " of " << fDstTree->GetEntries() << " in the file" << endl;
        ReadEvent(iEvent, "kalmanTracks");

        for (Int_t iPart = 0; iPart < fFemtoContainerReco->GetEntriesFast(); iPart++) {
            Int_t trId1 = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(iPart))->GetTrackID();

            for (Int_t jPart = iPart + 1; jPart < fFemtoContainerReco->GetEntriesFast(); jPart++) {
                Int_t trId2 = ((MpdFemtoContainer*) fFemtoContainerReco->UncheckedAt(jPart))->GetTrackID();

		// Nantes, 07 July

                cout << trId1 << " -- " << trId2 << endl;
                Float_t q = fCuts->Quality(trId1, trId2, iEvent);
//                Float_t s = fCuts->Sharing(trId1, trId2);
//
//                fHisto->GetQuality()->Fill(q);
//                if (fCuts->GetZeroSharing()) {
//                    fHisto->GetSharing()->Fill(s);
//                    fHisto->GetQualityVsSharing()->Fill(q, s);
//                } else {
//                    if (s != 0.) {
//                        fHisto->GetSharing()->Fill(s);
//                        fHisto->GetQualityVsSharing()->Fill(q, s);
//                    }
//                }
            }
        }
    }
}

Bool_t MpdFemto::Preselection(MpdTpcKalmanTrack* reco, MpdMCTrack* mc) {
    Float_t eta = -Log(Tan(reco->Theta() / 2.));

    if (mc->GetPdgCode() != fPDG)
        return kFALSE;

    else if (eta < fEtaCutLow || eta > fEtaCutUp)
        return kFALSE;

    else if (Abs(reco->Pt()) < fPtCutLow || Abs(reco->Pt()) > fPtCutUp)
        return kFALSE;

    else if (reco->GetNofHits() < fMinNoHits - 1)
        return kFALSE;

    else
        return kTRUE;
}
