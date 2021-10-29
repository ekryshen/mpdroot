/* 
 * A macro that illustrates how to work with the MpdFemto package
 * using different types of cuts (event, track, pair) and involving multiplicity analysis
 * DO NOT FORGET TO REDEINE A MULTIPLICITY MAP AS REQUIERED !!! 
 * It allows to construct 1d and 3d weighted correlation functions
 * IMPORTANT POINT: Input file passed as first argument must be written into the McDst-format
 * ($VMCWORKDIR/mpddst/mcDst)
 * ($VMCWORKDIR/macro/physical_analysis/createMcDstFromVHLLE.C)
 */

#include <Rtypes.h>
#include <TChain.h>
#include <TFile.h>
#include <cfloat>

using namespace std;

R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "macro/mpd/mpdloadlibs.C"

double M_ELECTRON = 0.000511;
double M_PION_PLUS = 0.139571;
double M_KAON_PLUS = 0.493677;
double M_PROTON = 0.938272;

// Prototypes of external functions to be used ...
void fillLikeSignAnalysis(MpdFemtoVertexMultAnalysis*, MpdFemtoBasicEventCut*, MpdFemtoBasicTrackCut*,
        MpdFemtoBasicPairCut*, MpdFemtoModelQinvCorrFctn*, MpdFemtoModelBPLCMS3DCorrFctnKt*);
void fillEventCut(MpdFemtoBasicEventCut*, MpdFemtoFxtEventCutMonitor*, MpdFemtoFxtEventCutMonitor*);
void fillTrackCut(MpdFemtoBasicTrackCut*, Int_t, Int_t, Int_t, MpdFemtoFxtTrackCutMonitor*, MpdFemtoFxtTrackCutMonitor*);
void fillTrackPairCut(MpdFemtoBasicPairCut*, MpdFemtoFxtPairCutMonitor*, MpdFemtoFxtPairCutMonitor*);

//________________

void hbtMcDstAdvanced(TString inFileList = "",
        TString oFileName = "oTest.root", Int_t nevents = 25) {
    
    map <Int_t, pair <Int_t, Int_t>> multiplicityBins;
    // bin --> (multLow, multUp)
    multiplicityBins[0] = make_pair(100, 120);
    multiplicityBins[1] = make_pair(120, 140);
    multiplicityBins[2] = make_pair(140, 160);
    multiplicityBins[3] = make_pair(160, 180);
    multiplicityBins[4] = make_pair(180, 200);
    
    cout << "Initializing McDstReader..." << endl;
    McDstReader* mcDstReader = new McDstReader(inFileList.Data());
    mcDstReader->Init();
    mcDstReader->setStatus("*", 0);
    mcDstReader->setStatus("Event*", 1);
    mcDstReader->setStatus("Particle*", 1);
    cout << "McDstReader initialized" << endl;

    if (!mcDstReader->chain()) {
        cout << "[ERROR] No chain has been found in McDst. Terminating." << endl;
        return;
    }

    Long64_t events2read = mcDstReader->chain()->GetEntries();
    cout << "Number of events in chain: " << events2read << endl;

    // Correct the amount event to read
    if (events2read >= nevents)
        events2read = nevents;

    // Create MpdFemtoMaker
    cout << "Initializing MpdFemtoMaker" << endl;
    MpdFemtoMaker* hbtMaker = new MpdFemtoMaker("HBT", "title");

    // Create MpdFemtoManager manager
    MpdFemtoManager* hbtManager = hbtMaker->hbtManager();

    // Setup MpdFemtoReader
    MpdFemtoMcDstReader *mcHbtReader = new MpdFemtoMcDstReader(mcDstReader);
    mcHbtReader->setRotateEventPlane(kFALSE);
    mcHbtReader->setMagneticField(-5.); // in kilogaus = -0.5 Tesla

    // Add reader to the manager
    hbtManager->setEventReader(mcHbtReader);

    // Lednicky weight generator                                                                                                                                                                                                             
    MpdFemtoModelWeightGeneratorLednicky* hbtWeight = new MpdFemtoModelWeightGeneratorLednicky();
    hbtWeight->setPairType(MpdFemtoModelWeightGeneratorLednicky::PionPlusPionPlus());
    hbtWeight->setCoulOff();
    hbtWeight->setQuantumOn();
    hbtWeight->setStrongOff();
    hbtWeight->set3BodyOff();

    // Theoretical model manager                                                                                                                                                                                                             
    MpdFemtoModelManager* thModelManager = new MpdFemtoModelManager();
    thModelManager->setWeightGenerator(hbtWeight);
    thModelManager->createCopyHiddenInfo(kTRUE);

    const Int_t nChargeBins = 2; // 0 - positive, 1 - negative
    TString charges[nChargeBins] = {"Positive", "Negative"};  

    const Int_t nMultBins = multiplicityBins.size();

    // Multiplicity analysis ... 
    MpdFemtoVertexMultAnalysis*** hbtAnalysis = new MpdFemtoVertexMultAnalysis**[nChargeBins];

    // 3D correlation functions ...
    MpdFemtoModelBPLCMS3DCorrFctnKt*** bpCorrFctn = new MpdFemtoModelBPLCMS3DCorrFctnKt**[nChargeBins];

    // 1D correlation functions ...
    MpdFemtoModelQinvCorrFctn*** qInvCorrFctn = new MpdFemtoModelQinvCorrFctn**[nChargeBins];

    // Prepare pointers for event, track and pair cuts ..  
    MpdFemtoBasicEventCut*** hbtEventCut = new MpdFemtoBasicEventCut**[nChargeBins];
    MpdFemtoBasicTrackCut*** hbtTrackCut = new MpdFemtoBasicTrackCut**[nChargeBins];
    MpdFemtoBasicPairCut*** hbtPairCut = new MpdFemtoBasicPairCut**[nChargeBins];

    // Loop over charges and chosen multiplicity ranges ...

    for (Int_t iCharge = 0; iCharge < nChargeBins; iCharge++) {
        hbtAnalysis[iCharge] = new MpdFemtoVertexMultAnalysis*[nMultBins];
        bpCorrFctn[iCharge] = new MpdFemtoModelBPLCMS3DCorrFctnKt*[nMultBins];
        qInvCorrFctn[iCharge] = new MpdFemtoModelQinvCorrFctn*[nMultBins];
        hbtEventCut[iCharge] = new MpdFemtoBasicEventCut*[nMultBins];
        hbtTrackCut[iCharge] = new MpdFemtoBasicTrackCut*[nMultBins];
        hbtPairCut[iCharge] = new MpdFemtoBasicPairCut*[nMultBins];

        for (auto it : multiplicityBins) {
            Int_t iMult = it.first;
            hbtEventCut[iCharge][iMult] = new MpdFemtoBasicEventCut();
            hbtTrackCut[iCharge][iMult] = new MpdFemtoBasicTrackCut();
            hbtPairCut[iCharge][iMult] = new MpdFemtoBasicPairCut();

            Int_t mEvents2Mix = 5;
            // Define local charge
            Int_t mLocalCharge = (iCharge == 0) ? 1 : -1;

            // Create vertex mult analysis
            Int_t multLow = it.second.first;
            Int_t multUp = it.second.second;
            hbtAnalysis[iCharge][iMult] = new MpdFemtoVertexMultAnalysis(1, -200., 200., 1, multLow, multUp);
            hbtAnalysis[iCharge][iMult]->setNumEventsToMix(mEvents2Mix);

            // Fill event cut and monitors                
            fillEventCut(hbtEventCut[iCharge][iMult], nullptr, nullptr);

            // Fill track cut and monitors                
            fillTrackCut(hbtTrackCut[iCharge][iMult], mLocalCharge, MpdFemtoBasicTrackCut::HbtPID::Pion, 0, nullptr, nullptr);

            // Fill pair cut and monitors
            fillTrackPairCut(hbtPairCut[iCharge][iMult], nullptr, nullptr);

            // Define 1D-analysis ... 
            Int_t qInvNbins = 200;
            double qInvRange[2] = {0., 2.};
            qInvCorrFctn[iCharge][iMult] =
                    new MpdFemtoModelQinvCorrFctn(Form("qinv_%sPart_multRange_%d_%d", charges[iCharge].Data(), multLow, multUp),
                    qInvNbins, qInvRange[0], qInvRange[1]);
            qInvCorrFctn[iCharge][iMult]->connectToManager(thModelManager);

            // Define 3D-analysis ...
            Int_t kTBins2 = 9;
            double kTRange2[2] = {0.15, 1.05};
            Int_t bpNbins = 40;
            double bpRange[2] = {-0.2, 0.2};

            bpCorrFctn[iCharge][iMult] = new MpdFemtoModelBPLCMS3DCorrFctnKt(Form("bp_%sPart_multRange_%d_%d", charges[iCharge].Data(), multLow, multUp),
                    bpNbins, bpRange[0], bpRange[1],
                    kTBins2, kTRange2[0], kTRange2[1]);
            bpCorrFctn[iCharge][iMult]->connectToManager(thModelManager);

            // Fill like sign analysis
            fillLikeSignAnalysis(hbtAnalysis[iCharge][iMult],
                    hbtEventCut[iCharge][iMult],
                    hbtTrackCut[iCharge][iMult],
                    hbtPairCut[iCharge][iMult],
                    qInvCorrFctn[iCharge][iMult],
                    bpCorrFctn[iCharge][iMult]);

            // Add multiplicity analysis ...
            hbtManager->addAnalysis(hbtAnalysis[iCharge][iMult]);
        }
    }

    Int_t iret = hbtMaker->init();

    // Loop over events
    for (Int_t iEvent = 0; iEvent < events2read; iEvent++) {
        cout << "Working on event# " << iEvent << endl;

        hbtMaker->clear();
        iret = hbtMaker->make();

        if (iret != 0) {
            cout << "[WARNING] hbt_dau_exp returned status: " << iret << ". Processing will be finished" << endl;
            break;
        }
    }

    // Write histograms to the output file
    TFile* oFile = new TFile(oFileName.Data(), "recreate");

    // Write histograms to the output file
    for (Int_t iCharge = 0; iCharge < nChargeBins; iCharge++) {
        for (Int_t iMult = 0; iMult < nMultBins; iMult++) {
            qInvCorrFctn[iCharge][iMult]->writeOutHistos();
            bpCorrFctn[iCharge][iMult]->writeOutHistos();
        }
    }
    oFile->Close();
}

// =============================================================

void fillLikeSignAnalysis(MpdFemtoVertexMultAnalysis* analysis,
        MpdFemtoBasicEventCut* eventCut,
        MpdFemtoBasicTrackCut* trackCut,
        MpdFemtoBasicPairCut* pairCut,
        MpdFemtoModelQinvCorrFctn* qInvCorrFctn,
        MpdFemtoModelBPLCMS3DCorrFctnKt* bpCorrFctn) {

    if (!analysis) {
        cout << "[ERROR] fillLikeSignAnalysis -- No analysis was provided \n";
        return;
    }

    if (!eventCut) {
        cout << "[ERROR] fillLikeSignAnalysis -- No event cut was provided \n";
        return;
    }

    if (!trackCut) {
        cout << "[ERROR] fillLikeSignAnalysis -- No track cut was provided \n";
        return;
    }

    if (!pairCut) {
        cout << "[ERROR] fillLikeSignAnalysis -- No pair cut was provided \n";
        return;
    }

    // Talk or not to talk
    analysis->setVerboseMode(kFALSE);

    // Set event cut
    analysis->setEventCut(eventCut);

    // Set first and second particle cuts
    analysis->setFirstParticleCut(trackCut);
    analysis->setSecondParticleCut(trackCut);

    // Set pair cut
    analysis->setPairCut(pairCut);

    // Set qInv correlation function if exists
    if (qInvCorrFctn)
        analysis->addCorrFctn(qInvCorrFctn);

    // Set Bertsch-Pratt correlation function if exists
    if (bpCorrFctn)
        analysis->addCorrFctn(bpCorrFctn);
}

void fillEventCut(MpdFemtoBasicEventCut *eventCut,
        MpdFemtoFxtEventCutMonitor *moniPass,
        MpdFemtoFxtEventCutMonitor *moniFail) {
    if (!eventCut) {
        cout << "[ERROR] fillEventCut -- No event cut was provided \n";
        return;
    }

    // Cut values
    float mImpact[2] = {-100., 100.};
    float vtxXShift = 0.;
    float vtxYShift = 0.;
    float vpdVzDiff[2] = {-5., 5.};
    float r[2] = {0., 2.};
    float sphericity[2] = {-1e6, 1e6};
    Int_t bTofMatch[2] = {0, 5000};
    Int_t bTofMult[2] = {0, 5000};
    Bool_t verbose = kFALSE;

    eventCut->setCheckBadRun(kFALSE);
    eventCut->setEventMult(0, 2000);
    eventCut->setImpactParameter(-1.1, -0.9); // A very important point when working with mcDst input got from vHLLE+UrQMD with no impact parameter written to the tree (-1 by default)
    eventCut->setVpdVzDiff(vpdVzDiff[0], vpdVzDiff[1]);
    eventCut->setVertXShift(vtxXShift);
    eventCut->setVertYShift(vtxYShift);
    eventCut->setVertRPos(r[0], r[1]);
    eventCut->setSphericity(sphericity[0], sphericity[1]); // Should not matter
    eventCut->setBTofTrayMult(bTofMult[0], bTofMult[1]);
    eventCut->setBTofMatchMult(bTofMatch[0], bTofMatch[1]);
    eventCut->setCent9(-1., 15);
    eventCut->setVerbose(verbose);

    if (moniPass && moniFail) {
        eventCut->addCutMonitor(moniPass, moniFail);
    }
}

//_________________

void fillTrackCut(MpdFemtoBasicTrackCut *trackCut,
        Int_t charge, Int_t particleType, Int_t detSel,
        MpdFemtoFxtTrackCutMonitor *moniPass,
        MpdFemtoFxtTrackCutMonitor *moniFail) {

    if (!trackCut) {
        cout << "[ERROR] fillTrackCut -- No track cut was provided\n";
        return;
    }

    // Standard cuts and kinematics
    Bool_t isTheory = kTRUE;
    Int_t pdgId = (charge > 0) ? 211 : -211;
    Bool_t primaryTracks = kTRUE;
    Bool_t verbose = kFALSE;
    Int_t nHits[2] = {15, 90};
    float nHitsRatio = 0.51;
    float pT[2] = {0.15, 2.8};
    float p[2] = {0., 2.8};
    float eta[2] = {-1., 1.};
    float dca[2] = {-0.1, 3.};

    // TPC identification only cuts
    float tpcMom[2] = {0.15, 0.80};
    float nSigmaTpcElectron[2] = {-2., 2.};
    float nSigmaTpcPion[2] = {-2., 2.};
    float nSigmaTpcKaon[2] = {-2., 2.};
    float nSigmaTpcProton[2] = {-2., 2.};
    float nSigmaTpcOther[2] = {-2., 2.};

    // TOF identification only cuts
    float tofMom[2] = {0.15, 1.5};
    float massSqr[2] = {-0.05, 0.08}; // for pions
    float invBetaDiff[2] = {-0.015, 0.015};

    // TPC+TOF identification
    float nSigmaTnTElectron[2] = {-3., 3.};
    float nSigmaTnTPion[2] = {-3., 3.};
    float nSigmaTnTKaon[2] = {-3., 3.};
    float nSigmaTnTProton[2] = {-3., 3.};

    // ToT with threshold
    float mPThresh = 0.45;

    // Set track cut parameters
    trackCut->selectPrimary(primaryTracks);
    trackCut->setCharge(charge);
    trackCut->setNHits(nHits[0], nHits[1]);
    trackCut->setAntiSplit(nHitsRatio);
    trackCut->setPt(pT[0], pT[1]);
    trackCut->setP(p[0], p[1]);
    trackCut->setEta(eta[0], eta[1]);
    trackCut->setDCA(dca[0], dca[1]);

    // Detector selection
    trackCut->setDetectorSelection(detSel);

    // TPC identification
    trackCut->setNSigmaElectron(nSigmaTpcElectron[0], nSigmaTpcElectron[1]);
    trackCut->setNSigmaPion(nSigmaTpcPion[0], nSigmaTpcPion[1]);
    trackCut->setNSigmaKaon(nSigmaTpcKaon[0], nSigmaTpcKaon[1]);
    trackCut->setNSigmaProton(nSigmaTpcProton[0], nSigmaTpcProton[1]);
    trackCut->setNSigmaOther(nSigmaTpcOther[0], nSigmaTpcOther[1]);
    trackCut->setTpcMom(tpcMom[0], tpcMom[1]);

    // TOF identification
    trackCut->setMassSqr(massSqr[0], massSqr[1]);
    trackCut->setInvBetaDiff(invBetaDiff[0], invBetaDiff[1]);
    trackCut->setTofMom(tofMom[0], tofMom[1]);

    // TPC+TOF identification
    trackCut->setTnTNSigmaElectron(nSigmaTnTElectron[0], nSigmaTnTElectron[1]);
    trackCut->setTnTNSigmaPion(nSigmaTnTPion[0], nSigmaTnTPion[1]);
    trackCut->setTnTNSigmaKaon(nSigmaTnTKaon[0], nSigmaTnTKaon[1]);
    trackCut->setTnTNSigmaProton(nSigmaTnTProton[0], nSigmaTnTProton[1]);

    // Set which particle type to select
    if (particleType == 1) {
        trackCut->setHbtPid(MpdFemtoBasicTrackCut::HbtPID::Electron);
        trackCut->setMass(M_ELECTRON);
    } else if (particleType == 2) {
        trackCut->setHbtPid(MpdFemtoBasicTrackCut::HbtPID::Pion);
        trackCut->setMass(M_PION_PLUS);
    } else if (particleType == 3) {
        trackCut->setHbtPid(MpdFemtoBasicTrackCut::HbtPID::Kaon);
        trackCut->setMass(M_KAON_PLUS);
    } else if (particleType == 4) {
        trackCut->setHbtPid(MpdFemtoBasicTrackCut::HbtPID::Proton);
        trackCut->setMass(M_PROTON);
    } else {
        cout << "[WARNING] fillTrackCut -- Unknown particle type: " << particleType <<
                "\nSwitching to pion PID" << endl;
        trackCut->setHbtPid(MpdFemtoBasicTrackCut::HbtPID::Pion);
        trackCut->setMass(M_PION_PLUS);
    }

    // Set momentum threshold for detector selection == 4
    trackCut->setPthresh(mPThresh);

    trackCut->setIsTheory(isTheory);
    trackCut->setPdgId(pdgId);

    // Print cut information for each track
    trackCut->setVerboseMode(verbose);

    // Add track cut monitors
    if (moniPass && moniFail) {
        trackCut->addCutMonitor(moniPass, moniFail);
    }
}

void fillTrackPairCut(MpdFemtoBasicPairCut *pairCut,
        MpdFemtoFxtPairCutMonitor *moniPass,
        MpdFemtoFxtPairCutMonitor *moniFail) {

    if (!pairCut) {
        cout << "[ERROR] fillTrackPairCut -- No pair cut was provided\n";
        return;
    }

    // Cut values
    float quality[2] = {-0.5, 0.6}; // antisplitting
    float kT[2] = {0.15, 1.05};
    float pT[2] = {0., 1e5};
    float openingAngle[2] = {-1e6, 1e6};
    float entranceSeparation[2] = {0., 1e6};
    float y[2] = {-1e5, 1e5};
    float eta[2] = {-1.5, 1.5};
    float qInv[2] = {-1e5, 1e5};
    float mInv[2] = {0., 1e5};
    float emissionAngle[2] = {-1e6, 1e6};
    float fmr[2] = {-1.1, 0.1};
    float closestRowAtDCA[2] = {-1e6, 1e6};
    float weightedAveSep[2] = {-1e6, 1e6};
    float aveSep[2] = {10., 1000.};
    float rValue = -1e6;
    float dPhiStarMin[2] = {-1e6, 1e6};
    Bool_t verbose = kFALSE;

    // Set pair cut values
    pairCut->setQuality(quality[0], quality[1]);
    pairCut->setKt(kT[0], kT[1]);
    pairCut->setPt(pT[0], pT[1]);
    pairCut->setOpeningAngle(openingAngle[0], openingAngle[1]);
    pairCut->setEntranceSeparation(entranceSeparation[0], entranceSeparation[1]);
    pairCut->setRapidity(y[0], y[1]);
    pairCut->setEta(eta[0], eta[1]);
    pairCut->setQinv(qInv[0], qInv[1]);
    pairCut->setMinv(mInv[0], mInv[1]);
    pairCut->setAngleToPrimaryVertex(emissionAngle[0], emissionAngle[1]);
    pairCut->setFracOfMergedRow(fmr[0], fmr[1]);
    pairCut->setClosestRowAtDCA(closestRowAtDCA[0], closestRowAtDCA[1]);
    pairCut->setWeightedAvSep(weightedAveSep[0], weightedAveSep[1]);
    pairCut->setAverageSeparation(aveSep[0], aveSep[1]);
    pairCut->setRValue(rValue);
    pairCut->setDPhiStarMin(dPhiStarMin[0], dPhiStarMin[1]);
    pairCut->setVerboseMode(verbose);

    // Add pair cut monitors
    if (moniPass && moniFail) {
        pairCut->addCutMonitor(moniPass, moniFail);
    }
}
