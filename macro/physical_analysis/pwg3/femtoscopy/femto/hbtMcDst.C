/* 
 * A macro that illustrates how to work with the MpdFemto package
 * It allows to construct 1d and 3d weighted correlation functions
 * IMPORTANT POINT: Input file passed as first argument must be written into the McDst-format
 * ($VMCWORKDIR/mpddst/mcDst)
 * ($VMCWORKDIR/macro/physical_analysis/createMcDstFromVHLLE.C)
 */

#include <Rtypes.h>
#include <TChain.h>
#include <TFile.h>

using namespace std;

R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "macro/mpd/mpdloadlibs.C"

void hbtMcDst(TString fileName = "",
        TString oFileName = "oTest.root", Int_t nevents = 100) {

    cout << " fileName: " << fileName << endl;
    cout << " nevents:  " << nevents << endl;

    //
    // Monte Carlo input setup
    //
    McDstReader *mcDstReader = new McDstReader(fileName);
    mcDstReader->Init();
    mcDstReader->setStatus("*", 0);
    mcDstReader->setStatus("Event", 1);
    mcDstReader->setStatus("Particle", 1);

    if (!mcDstReader->chain())
        cout << "No chain has been found." << endl;

    Long64_t events2read = mcDstReader->chain()->GetEntries();
    cout << "Number of events in chain: " << events2read << endl;

    // Correct the amount event to read
    if (events2read >= nevents)
        events2read = nevents;

    //
    // MpdFemtoMaker setup
    //

    // MpdFemtoMaker
    MpdFemtoMaker* hbtMaker = new MpdFemtoMaker("HBT", "title");
    cout << "MpdFemtoMaker instantiated" << endl;

    // Create manager
    MpdFemtoManager* hbtManager = hbtMaker->hbtManager();

    // Setup MpdFemtoReader
    MpdFemtoMcDstReader *mcHbtReader = new MpdFemtoMcDstReader(mcDstReader);
    mcHbtReader->setRotateEventPlane(false);
    mcHbtReader->setMagneticField(-5.); // in kilogaus = -0.5 Tesla

    // Add reader to the manager
    hbtManager->setEventReader(mcHbtReader);

    // Create and set analysis
    MpdFemtoAnalysis *hbtAnalysis = new MpdFemtoAnalysis();
    hbtAnalysis->setVerboseMode(false);

    // Create and set event cut
    MpdFemtoBasicEventCut *eventCut = new MpdFemtoBasicEventCut();
    // A very important point when working with mcDst input got from vHLLE+UrQMD with no impact parameter written to the tree (-1 by default)
    eventCut->setImpactParameter(-1.1, -0.9);
    hbtAnalysis->setEventCut(eventCut);

    // Create and set track cut
    MpdFemtoBasicTrackCut *trackCut = new MpdFemtoBasicTrackCut();
    trackCut->setIsTheory(true);
    trackCut->setPdgId(211);
    trackCut->setEta(-1., 1.);
    trackCut->setPt(0.15, 1.55);
    trackCut->setMass(M_PION_PLUS); // from MpdFemtoMaker/phys_constants.h

    // Create and set cut monitor
    MpdFemtoFxtTrackCutMonitor *trackMoniPass = new MpdFemtoFxtTrackCutMonitor("_moniPass_", M_PION_PLUS);
    MpdFemtoFxtTrackCutMonitor *trackMoniFail = new MpdFemtoFxtTrackCutMonitor("_moniFail_", M_PION_PLUS);
    trackCut->addCutMonitor(trackMoniPass, trackMoniFail);

    // Set particle cuts
    hbtAnalysis->setFirstParticleCut(trackCut);
    hbtAnalysis->setSecondParticleCut(trackCut);

    // Create and set pair cut
    MpdFemtoBasicPairCut *pairCut = new MpdFemtoBasicPairCut();

    // Create pair cut monitors
    MpdFemtoFxtPairCutMonitor *pairMoniPass = new MpdFemtoFxtPairCutMonitor("_moniPass_");
    MpdFemtoFxtPairCutMonitor *pairMoniFail = new MpdFemtoFxtPairCutMonitor("_moniFail_");
    pairCut->addCutMonitor(pairMoniPass, pairMoniFail);

    // Add pair cut to the analysis
    hbtAnalysis->setPairCut(pairCut);

    // Set how many events to mix
    hbtAnalysis->setNumEventsToMix(5);

    // Lednicky weight generator
    MpdFemtoModelWeightGeneratorLednicky *hbtWeight = new MpdFemtoModelWeightGeneratorLednicky();
    hbtWeight->setPairType(MpdFemtoModelWeightGeneratorLednicky::PionPlusPionPlus());
    hbtWeight->setCoulOff();
    hbtWeight->setQuantumOn();
    hbtWeight->setStrongOff();
    hbtWeight->set3BodyOff();

    // Theoretical model manager
    MpdFemtoModelManager *thModelManager = new MpdFemtoModelManager();
    thModelManager->setWeightGenerator(hbtWeight);
    thModelManager->createCopyHiddenInfo(true);

    // Create 1D correlation function integrated over kT
    /// Ordered list of parameters to be passed to the constructor ...
    /// \param title  Name of the histogram
    /// \param nBins  Number of bins
    /// \param  qInvLow   Minimum value of the qInv
    /// \param  qInvHi    Minimum value of the qInv

    MpdFemtoModelQinvCorrFctn *thCorrFctn = new MpdFemtoModelQinvCorrFctn("hTheorQinv", 40, 0., 0.4);
    thCorrFctn->connectToManager(thModelManager);

    // Create 3D correlation function integrated with kT binning 
    /// Ordered list of parameters to be passed to the constructor ... 
    /// \param title  Name of the histogram
    /// \param nBins  Number of bins (will be used for out, side and long projections)
    /// \param qLo    Minimum value of the q
    /// \param qHi    Maximum value of the q
    /// \param ktBins Number of kT bins used in the analysis
    /// \param ktLo   Minimum value of kT
    /// \param ktHi   Maximum value of kT
    MpdFemtoModelBPLCMS3DCorrFctnKt *thCorrFctn3D = new MpdFemtoModelBPLCMS3DCorrFctnKt("hTheorBPLCMS", 80, -0.4, 0.4, 9, 0.05, 0.95);
    thCorrFctn3D->connectToManager(thModelManager);

    // Add correlation function to the analysis
    hbtAnalysis->addCorrFctn(thCorrFctn);
    hbtAnalysis->addCorrFctn(thCorrFctn3D);

    // Add analysis to the manager
    hbtManager->addAnalysis(hbtAnalysis);

    // Initialize MpdFemtoMaker
    int iret = hbtMaker->init();
    // Loop over events
    for (int iEvent = 0; iEvent < events2read; iEvent++) {

        hbtMaker->clear();
        cout << "hbtThQinv -- Working on event # " << iEvent << endl;
        iret = hbtMaker->make();

        if (iret != 0) {
            cout << "[WARNING] MpdFemtoMaker returned status: " << iret
                    << ". Processing will be finished" << endl;
            break;
        }
    }

    // Maker finish flag
    iret = hbtMaker->Finish();

    // Write histograms to the output file
    TFile *oFile = new TFile(oFileName.Data(), "recreate");
    trackMoniPass->writeOutHistos();
    trackMoniFail->writeOutHistos();
    pairMoniPass->writeOutHistos();
    pairMoniFail->writeOutHistos();
    thCorrFctn->writeOutHistos();
    thCorrFctn3D->writeOutHistos();
    oFile->Close();
}
