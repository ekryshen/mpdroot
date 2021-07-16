#include "MpdKfParticleFinder.h"

MpdKfParticleFinder::MpdKfParticleFinder(TString dstName) :
fTpcTracks(nullptr),
fMcTracks(nullptr),
fMiniTracks(nullptr),
fMiniCovMatrices(nullptr),
fPrimaryVertices(nullptr),
fTopoReconstructor(nullptr),
fPerformance(nullptr),
fSecondaries(nullptr),
fKFVertices(nullptr),
fMCHeader(nullptr),
fMCEvent(nullptr),
fHeader(nullptr),
isUseCuts(kFALSE),
isMini(kFALSE),
isPerformance(kFALSE),
fCuts(nullptr), fPdgDB(nullptr), fInput("") {
    fVerbose = 0;

    // Setting default PID hypo to be used ...
    fPidHypo.first = 2212;
    fPidHypo.second = -211;

    // Initialize FairRoot manager
    ioman = FairRootManager::Instance();

    FairFileSource* fFileSource = nullptr;
    FairSource* fSource = nullptr;
    MpdMiniDstFileSource* fMiniSource = nullptr;

    // It can be a single file or a list of files (*.lis or *.list) to process ...
    if (dstName.Contains(".lis") || dstName.Contains(".list")) {

        // Reading passed list with inputs ...
        ReadList(dstName);

        if (isMini) {
            fMiniSource = new MpdMiniDstFileSource(*fInputs.begin());
            for (auto it = next(fInputs.begin(), 1); it != fInputs.end(); it++)
                fMiniSource->AddFile(*it);

            fSource = (FairSource*) fMiniSource;
        } else {
            fFileSource = new FairFileSource(*fInputs.begin());

            for (auto it = next(fInputs.begin(), 1); it != fInputs.end(); it++)
                fFileSource->AddFile(*it);
        }
    } else if (dstName.Contains(".root")) {
        fInput = dstName;

        if (dstName.Contains("MiniDst")) {
            fSource = new MpdMiniDstFileSource(fInput);
            isMini = kTRUE;
        } else
            fFileSource = new FairFileSource(fInput);
    } else
        Fatal("MpdKfParticleFinder::MpdKfParticleFinder(TString dstName)", "Provided input seems to be incorrect!");

    // Setting already specified file source ...
    if (fFileSource)
        FairRunAna::Instance()->SetSource(fFileSource);
    else if (fSource)
        FairRunAna::Instance()->SetSource(fSource);
}

void MpdKfParticleFinder::ReadList(TString dstName) {
    string line;
    ifstream f(dstName.Data(), ios::in);

    TString dst = "";
    Int_t id = 0;

    const Int_t nMaxFiles = 5000;
    Int_t iCurrentFile = 0;

    while (!f.eof()) {
        f >> id >> dst;
        fInputs.push_back(dst);
        fIds[dst] = id;

        getline(f, line);
        iCurrentFile++;

        if (iCurrentFile > nMaxFiles)
            Fatal("MpdKfParticleFinder::ReadList()", "Exceeded number of files or something wrong with passed list!");
    }

    for (auto it = fInputs.begin(); it != fInputs.end(); it++) {
        dst = *it;

        if (!dst.Contains(".root"))
            fInputs.erase(it--);
    }

    auto it = fInputs.begin();
    isMini = ((*it).Contains("MiniDst") ? kTRUE : kFALSE);
}

InitStatus MpdKfParticleFinder::Init() {
    if (isUseCuts)
        fCuts = new TrackCuts();

    // Getting necessary branches to be used ...
    // It depends on a dst filename passed as argument

    if (isMini) {
        fMiniTracks = (TClonesArray*) ioman->GetObject("Track");
        fMiniCovMatrices = (TClonesArray*) ioman->GetObject("TrackCovMatrix");
        fMCEvent = (TClonesArray*) ioman->GetObject("McEvent");
    } else {
        // Tpc tracks ...
        fTpcTracks = (TClonesArray*) ioman->GetObject("TpcKalmanTrack");

        // Monte Carlo tracks ...
        fMcTracks = (TClonesArray*) ioman->GetObject("MCTrack");

        // Reconstructed PV's ...
        fPrimaryVertices = (TClonesArray*) ioman->GetObject("Vertex");

        // MC event header ...
        fMCHeader = (FairMCEventHeader*) ioman->GetObject("MCEventHeader.");
    }

    fHeader = (FairEventHeader*) ioman->GetObject("EventHeader.");

    // Register output branch ...
    fSecondaries = new TClonesArray("KFParticle");
    ioman->Register("Secondaries", "Secondaries_", fSecondaries, kTRUE);

    fKFVertices = new TClonesArray("KFVertex");
    ioman->Register("Vertex", "Vertex_", fKFVertices, kTRUE);

    // Calling topoReconstructor ...
    fTopoReconstructor = new KFParticleTopoReconstructor();
    fTopoReconstructor->SetField(5.); // FIXME !!! (to be read automatically)

    // Setting a decay to be reconstructed ...
    map <Int_t, Bool_t> decayList;

    Int_t decayPdg = -1;
    if (fPidHypo.first == 2212 && fPidHypo.second == -211)
        decayPdg = 3122; // \Lambda^{0} -> p + \pion^{-}
    else if (fPidHypo.first == 211 && fPidHypo.second == -211)
        decayPdg = 310; // \Kaon^{0}_{s} -> \pion^{+} + \pion^{-}

    if (decayPdg != -1)
        decayList[decayPdg] = kTRUE;

    fTopoReconstructor->GetKFParticleFinder()->SetReconstructionList(decayList);
    map <Int_t, Bool_t> decayListToCheck = fTopoReconstructor->GetKFParticleFinder()->GetReconstructionList();

    if (decayListToCheck.size() != 1)
        Fatal("MpdKfParticleFinder::Init()", "Only one decay in the list is supported !!!"); // FIXME !!! 

    for (auto pdg : decayListToCheck)
        fDecPdgs.push_back(pdg.first);

    // Let's see a list of decays to be reconstructed (if necessary) ...
    if (fVerbose)
        for (auto pdgDecay : decayListToCheck)
            cout << "PDG# " << pdgDecay.first << " isActive# " << pdgDecay.second << endl;

    // Right now KFPerformance does not support miniDst!
    if (isMini) 
        isPerformance = kFALSE;
    
    if (isPerformance) {
        // Involving performance tools ...
        fPerformance = new KFTopoPerformance();
        // fPerformance->DoNotStoreMCHistograms();

        TString outPerf = (TString) FairRunAna::Instance()->GetOutputFile()->GetName();
        outPerf = "KFPerformance_" + outPerf;
        fOutFile = new TFile(outPerf.Data(), "recreate");
        fPerformance->CreateHistos("KFHistos", (TDirectory*) fOutFile, decayListToCheck);

        fPdgDB = new TDatabasePDG();
    }

    return kSUCCESS;
}

void MpdKfParticleFinder::Exec(Option_t* option) {
    if (evCounter % 1000 == 0)
        cout << "Processing Event# " << evCounter << endl;

    if (fInputs.size()) {
        auto it = fIds.find((TString) ioman->GetRootFile()->GetName());
        MpdMiniMcEvent* mcEv = (MpdMiniMcEvent*) fMCEvent->UncheckedAt(0);

        if (isMini)
            fHeader->SetMCEntryNumber(mcEv->eventId());

        if (it != fIds.end())
            fHeader->SetInputFileId(it->second);
        else
            Fatal("MpdKfParticleFinder::Exec", "Input file is not associated with its id!");
    }

    fTopoReconstructor->Clear();
    fSecondaries->Delete();
    fKFVertices->Delete();

    ProcessDst();

    evCounter++;
}

void MpdKfParticleFinder::Finish() {
    if (isPerformance) {
        // Performance tools ...
        fPerformance->GetHistosDirectory()->Write();
        fOutFile->Close();
    }

    // PDG DB ...
    if (fPdgDB)
        delete fPdgDB;
}

void MpdKfParticleFinder::Mini2Kalman(TClonesArray* arr) {
    for (Int_t iTrack = 0; iTrack < fMiniTracks->GetEntriesFast(); iTrack++) {
        MpdMiniTrack* mini = (MpdMiniTrack*) fMiniTracks->UncheckedAt(iTrack);
        MpdMiniTrackCovMatrix* cMatrix = (MpdMiniTrackCovMatrix*) fMiniCovMatrices->UncheckedAt(iTrack);

        MpdTpcKalmanTrack* kalTrack = new ((*arr)[arr->GetEntriesFast()]) MpdTpcKalmanTrack();
        kalTrack->SetPos(cMatrix->R0());
        kalTrack->SetPosNew(cMatrix->R());
        kalTrack->SetNofHits(mini->nHits());

        TMatrixD stateVector = cMatrix->stateVector();
        kalTrack->SetParam(stateVector);

        TMatrixDSym covMatrix = cMatrix->covarianceMatrix();
        kalTrack->SetCovariance(covMatrix);
    }
}

void MpdKfParticleFinder::ProcessDst() {
    vector <ExtendedMpdTpcKalmanTrack> tpcTracks;

    TClonesArray* tpcKalmanTracks = nullptr;
    if (isMini) {
        tpcKalmanTracks = new TClonesArray("MpdTpcKalmanTrack");

        Mini2Kalman(tpcKalmanTracks);
    }

    // 1. Loop over found TPC tracks ...
    TClonesArray* inTrackArray = isMini ? tpcKalmanTracks : fTpcTracks;

    for (Int_t iTpcTrack = 0; iTpcTrack < inTrackArray->GetEntriesFast(); iTpcTrack++) {
        MpdTpcKalmanTrack* tr = (MpdTpcKalmanTrack*) inTrackArray->UncheckedAt(iTpcTrack);
        ExtendedMpdTpcKalmanTrack track = *(ExtendedMpdTpcKalmanTrack*) tr;

        if (isUseCuts) {
            // Skipping tracks not satisfying to the established cuts ...
            Double_t eta = Log(Tan(track.Theta() / 2.));
            if (track.Pt() < fCuts->GetPtMin() || track.Pt() > fCuts->GetPtMax() || Abs(eta) > fCuts->GetAbsEta())
                continue;
        }

        track.fIdx = iTpcTrack;
        tpcTracks.push_back(track);
    }

    if (fVerbose)
        cout << "Number of found TPC tracks# " << tpcTracks.size() << endl;

    if (tpcTracks.size() == 0)
        return;

    // 2. Getting and adding information on primary vertex ...
    vector <Int_t> pvTrackIds;

    if (isMini)
        for (Int_t iMiniTrack = 0; iMiniTrack < fMiniTracks->GetEntriesFast(); iMiniTrack++) {
            MpdMiniTrack* mini = (MpdMiniTrack*) fMiniTracks->UncheckedAt(iMiniTrack);

            if (mini->isPrimary())
                pvTrackIds.push_back(iMiniTrack);
        } 
    else {
        MpdVertex* recoVp = (MpdVertex*) fPrimaryVertices->UncheckedAt(0);
        Int_t* indcs = recoVp->GetIndices()->GetArray();
        Int_t size = recoVp->GetIndices()->GetSize();
        pvTrackIds.assign(indcs, indcs + size);
    }

    // 3. Filling KFPTrackVector (xyz, PxPyPz, Q, PDG) with data ...
    KFPTrackVector trVector;
    FillKFPTrackVector(tpcTracks, &trVector, pvTrackIds);

    // 4. Getting KFPTracks from KFPTrackVector ...
    vector <KFParticle> kfParticles;
    vector <Int_t>* pdgHypos = new vector <Int_t>;
    for (Int_t iTrack = 0; iTrack < trVector.Size(); iTrack++) {
        KFPTrack track;
        trVector.GetTrack(track, iTrack);
        Int_t q = track.Charge();
        Int_t hypo = (q > 0) ? fPidHypo.first : fPidHypo.second;

        KFParticle particle(track, hypo);
        particle.SetId(track.Id());

        kfParticles.push_back(particle);
        pdgHypos->push_back(hypo);
    }

    // 5. Initializing topoReconstructor calling one of possible user constructors ...
    fTopoReconstructor->Init(kfParticles, pdgHypos);

    // 6. Reconstruct primary vertex ...
    fTopoReconstructor->ReconstructPrimVertex();

    // 7. Reconstruct secondary particles and select best candidates ...
    fTopoReconstructor->SortTracks();
    fTopoReconstructor->ReconstructParticles();
    fTopoReconstructor->SelectParticleCandidates();

    // 8. Writing found candidates to output file ...
    for (KFParticle recoSecondary : fTopoReconstructor->GetParticles())
        new ((*fSecondaries) [fSecondaries->GetEntriesFast()]) KFParticle(recoSecondary);

    delete pdgHypos;

    // 9. Getting reconstructed primary vertex to write it to output ...
    KFVertex vertex = fTopoReconstructor->GetPrimKFVertex();
    new ((*fKFVertices) [fKFVertices->GetEntriesFast()]) KFVertex(vertex);

    if (isPerformance && !isMini) {
        // Involving KFParticlePerformance ...
        fPerformance->SetTopoReconstructor(fTopoReconstructor);

        // Filling KFMCTrack for current event ....
        vector <KFMCTrack> mcTracksInEvent = FillKFMCTrack();

        // Setting MC tracks ...
        fPerformance->SetMCTracks(mcTracksInEvent);

        // Checking MC tracks ...
        fPerformance->CheckMCTracks();

        // Doing KFParticle --> KFMCParticle (nDaughters = 1)
        vector <Int_t> kf2kfmc;
        kf2kfmc.resize(fTpcTracks->GetEntriesFast());

        for (Int_t iRP = 0; iRP < fTopoReconstructor->GetParticles().size(); iRP++) {
            const KFParticle& rPart = fTopoReconstructor->GetParticles()[iRP];

            if (rPart.NDaughters() != 1)
                continue;

            MpdKalmanTrack* tr = (MpdKalmanTrack*) fTpcTracks->UncheckedAt(rPart.DaughterIds().at(0));
            kf2kfmc.at(rPart.DaughterIds().at(0)) = tr->GetTrackID();
        }

        fPerformance->SetTrackMatch(kf2kfmc);

        // Matching mc and reco tracks ...
        fPerformance->MatchTracks();

        // Filling a set of QA histos ... 
        fPerformance->FillHistos();
    }

    if (tpcKalmanTracks)
        delete tpcKalmanTracks;
}

vector <KFMCTrack> MpdKfParticleFinder::FillKFMCTrack() {
    vector <KFMCTrack> tracks;
    tracks.clear();

    // Getting Monte Carlo ids of reconstructed tracks ...
    vector <Int_t> recoMcIdx;
    recoMcIdx.clear();

    for (Int_t iMc = 0; iMc < fMcTracks->GetEntriesFast(); iMc++) {
        MpdMCTrack* mcTrack = (MpdMCTrack*) fMcTracks->UncheckedAt(iMc);

        KFMCTrack track;

        track.SetX(mcTrack->GetStartX());
        track.SetY(mcTrack->GetStartY());
        track.SetZ(mcTrack->GetStartZ());

        track.SetPx(mcTrack->GetPx());
        track.SetPy(mcTrack->GetPy());
        track.SetPz(mcTrack->GetPz());

        track.SetPDG(mcTrack->GetPdgCode());
        track.SetMotherId(mcTrack->GetMotherId());

        TParticlePDG* pdgParticle = fPdgDB->GetParticle(mcTrack->GetPdgCode());

        Int_t q = (pdgParticle) ? pdgParticle->Charge() / 3 : 0;

        track.SetQP(q / mcTrack->GetP());
        track.SetNMCPoints(mcTrack->GetNPoints(kTPC));

        // Trying to know whether it was reconstructed ...
        auto it = find(recoMcIdx.begin(), recoMcIdx.end(), iMc);

        if (it != recoMcIdx.end())
            track.SetReconstructed();
        else
            track.SetNotReconstructed();

        // Is it in the TPC acceptance? ...        
        if (!isInAcceptance(mcTrack))
            track.SetOutOfDetector();

        tracks.push_back(track);
    }

    return tracks;
}

Bool_t MpdKfParticleFinder::isInAcceptance(MpdMCTrack* track) {
    Double_t eta = TMath::ATanH(track->GetPz() / track->GetP());
    Double_t pt = track->GetPt();

    if (TMath::Abs(eta) < 1.2 && pt > .15 && pt < 2.) // FIXME !!!
        return kTRUE;
    else
        return kFALSE;
}

void MpdKfParticleFinder::tpcInnerShell(MpdTpcKalmanTrack in, MpdTpcKalmanTrack& out) {
    // Let's get covariance matrix and track params. approximately corresponding to the inner TPC shell ...
    // We need to have a copy of considering track to be able to fill fParamNew
    // to be instantiated when doing the copy ... 
    MpdTpcKalmanTrack tmp(in);
    tmp.SetDirection(MpdKalmanTrack::kOutward);

    /// Propagate TPC track to 27 cm inward (fPos == 27)
    MpdKalmanHit hitTmp;
    hitTmp.SetType(MpdKalmanHit::kFixedR);
    hitTmp.SetPos(tmp.GetPos()); // fPos ~= 27 cm

    tmp.SetParamNew(*tmp.GetParam());
    tmp.SetPos(tmp.GetPosNew());
    tmp.ReSetWeight();
    TMatrixDSym w = *tmp.GetWeight(); // save current weight matrix

    Bool_t ok = MpdKalmanFilter::Instance()->PropagateToHit(&tmp, &hitTmp, kFALSE, kFALSE);

    if (!ok) {
        out = tmp;
        return;
    }
    // tmp.SetWeight(w); // restore original weight matrix (near TPC inner shell)

    // После них GetCovariance() и GetParamNew() дадут ков. матр. и параметры в точке GetPosNew(), которая около R=27 см (на внутренней оболочке TPC). 
    TMatrixDSym cov(*tmp.Weight2Cov());
    TMatrixD param(*tmp.GetParamNew());

    // Resetting params and cov. matrix by the propagated ones ...
    tmp.SetParam(param);
    tmp.SetCovariance(cov);

    out = tmp;
}

void MpdKfParticleFinder::lastHit(MpdTpcKalmanTrack in, MpdTpcKalmanTrack& out) {

    MpdTpcKalmanTrack tmp(in);

    // Resetting params and cov. matrix to those ones  at the last hit ...
    tmp.SetParam(*in.GetParamAtHit());

    TMatrixDSym* weightAtLastHit = in.GetWeightAtHit();
    TMatrixDSym& covarAtLastHit = weightAtLastHit->InvertFast();

    tmp.SetCovariance(covarAtLastHit);

    out = tmp;
}

void MpdKfParticleFinder::dcaToBeamline(MpdTpcKalmanTrack in, MpdTpcKalmanTrack& out) {

    MpdTpcKalmanTrack tmp(in);

    // Resetting params and cov. matrix to those ones at dca2beamline ...
    tmp.SetParam(*in.GetParam());
    tmp.SetCovariance(*in.GetCovariance());

    out = tmp;
}

void MpdKfParticleFinder::FillKFPTrackVector(vector <ExtendedMpdTpcKalmanTrack> tracksFromTpc, KFPTrackVector* tracks, vector <Int_t> trIdxPv) {
    tracks->Resize(tracksFromTpc.size());

    for (Int_t iTrack = 0; iTrack < tracksFromTpc.size(); iTrack++) {
        ExtendedMpdTpcKalmanTrack in = tracksFromTpc[iTrack];

        MpdTpcKalmanTrack out;

        if (fVerbose) {
            cout << "Track params. at 2D-DCA to beamline: " << endl;
            in.GetParam()->Print();
            in.GetCovariance()->Print();
        }

        // Parameters are propagated to the inner TPC shell at R of order of 27 cm
        tpcInnerShell(in, out);
        // lastHit(in, out);
        // dcaToBeamline(in, out);

        if (fVerbose) {
            cout << "Propagated parameters: " << endl;
            out.GetParam()->Print();
            out.GetCovariance()->Print();
        }

        // Getting propagated track coordinates / angles / momenta  ...
        Double_t phi = out.GetParam(0) / out.GetPosNew();
        Double_t X = out.GetPosNew() * Cos(phi);
        Double_t Y = out.GetPosNew() * Sin(phi);
        Double_t Z = out.GetZ();

        Double_t q = out.Charge();

        Double_t Px = out.Momentum3().X();
        Double_t Py = out.Momentum3().Y();
        Double_t Pz = out.Momentum3().Z();

        // Setting state vector ....
        vector <Double_t> parVectorToSet{X, Y, Z, Px, Py, Pz};
        for (Int_t iParam = 0; iParam < parVectorToSet.size(); iParam++)
            tracks->SetParameter(parVectorToSet[iParam], iParam, iTrack);

        tracks->SetQ(q, iTrack);
        tracks->SetId(in.fIdx, iTrack);

        // Setting PID hypo ...
        if (q > 0)
            tracks->SetPDG(fPidHypo.first, iTrack);
        else
            tracks->SetPDG(fPidHypo.second, iTrack);

        // Setting PVIndex = 0 to tracks considered as primaries ...
        auto it = find(trIdxPv.begin(), trIdxPv.end(), iTrack);
        if (it != trIdxPv.end())
            tracks->SetPVIndex(0, iTrack);
            // Probably, looks as a secondary ...
        else
            tracks->SetPVIndex(-1, iTrack);

        vector <Double_t> cov21 = DoErrorPropagationToXYZPxPyPz(out.GetCovariance(), out.GetParam(), &out);

        for (Int_t iCov = 0; iCov < cov21.size(); iCov++)
            tracks->SetCovariance(cov21.at(iCov), iCov, iTrack);
    }
}

vector <Double_t> MpdKfParticleFinder::DoErrorPropagationToXYZPxPyPz(TMatrixDSym* cov, TMatrixD * params, MpdTpcKalmanTrack* track) {

    // State track vector in MPD: 
    // 0) r\Phi_{T} 
    // 1) z 
    // 2) \Phi 
    // 3) \Lambda (dip angle) 
    // 4) -q / Pt

    // J : A --> B 
    // A : x, y, z, px, py, pz 
    // B : x, y, z, \Phi, \Lambda, -q / Pt

    // Getting Jacobian ...
    Float_t J[6][6];
    for (Int_t i = 0; i < 6; i++)
        for (Int_t j = 0; j < 6; j++)
            J[i][j] = 0;

    // Getting corresponding track parameters ...
    Int_t q = track->Charge();
    Double_t pt = track->Pt();

    Double_t Phi = (*params)(2, 0);
    Double_t Lambda = (*params)(3, 0);

    // Coordinate transformations ...
    J[0][0] = 1.; // dx / dx
    J[1][1] = 1.; // dy / dy
    J[2][2] = 1.; // dz / dz

    // Momentum transformations ...
    J[3][3] = -pt * Sin(Phi); // dPx / d\Phi
    J[3][5] = pt * pt * Cos(Phi) / q; // dPx / d(-q / Pt)

    J[4][3] = pt * Cos(Phi); // dPy / d\Phi
    J[4][5] = pt * pt * Sin(Phi) / q; // dPy / d(-q / Pt)

    J[5][4] = pt / (Cos(Lambda) * Cos(Lambda)); // dPz / dLambda
    J[5][5] = Tan(Lambda) * pt * pt / q; // dPz / d(-q / pt)

    // Extending track covariance matrix by one row and one column ...
    Float_t CovIn[6][6]; // triangular -> symmetric matrix

    CovIn[0][0] = 1e-4; // dx. start from nowhere

    for (Int_t i = 1; i < 6; i++) {
        CovIn[i][0] = 0;
        CovIn[0][i] = 0;
    }

    for (Int_t i = 1; i < 6; i++) {
        for (Int_t j = 1; j <= i; j++) {
            CovIn[i][j] = (*cov)(i - 1, j - 1);
            CovIn[j][i] = (*cov)(i - 1, j - 1);
        }
    }

    Float_t CovInJt[6][6]; // CovInJt = CovIn * J^t
    for (Int_t i = 0; i < 6; i++)
        for (Int_t j = 0; j < 6; j++) {
            CovInJt[i][j] = 0;
            for (Int_t k = 0; k < 6; k++)
                CovInJt[i][j] += CovIn[i][k] * J[j][k];
        }

    Float_t CovOut[6][6]; // CovOut = J * CovInJt
    for (Int_t i = 0; i < 6; i++)
        for (Int_t j = 0; j < 6; j++) {
            CovOut[i][j] = 0;
            for (Int_t k = 0; k < 6; k++)
                CovOut[i][j] += J[i][k] * CovInJt[k][j];
        }

    vector <Double_t> KFPCov; // symmetric matrix -> triangular

    for (Int_t i = 0; i < 6; i++)
        for (Int_t j = 0; j <= i; j++)
            KFPCov.push_back(CovOut[i][j]);

    return KFPCov;
}