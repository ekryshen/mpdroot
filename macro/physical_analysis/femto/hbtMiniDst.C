/* 
 * A macro that illustrates how to work with the MpdFemtoMaker package
 * It allows to read and analyse data in MiniDst format.
 * As an input file the macro accepts either fname.MiniDst.root file or
 * a file list fname.list (or fname.lis) that contains a list of fname.MiniDst.root files.
 * IMPORTANT POINT: The files that will be analysed MUST have .MiniDst.root extension.
 */

#include "Rtypes.h"
#include "TChain.h"
#include "TFile.h"
#include "TROOT.h"
#include "TObject.h"
#include "Riostream.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH3.h"

// C++ headers
#include <iostream>

R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "macro/mpd/mpdloadlibs.C"

//________________
void hbtMiniDst(const char* inFileList = "",
		const char* oFileName = "oHbtMiniDst.root") {

  std::cout << " Input file list: " << inFileList << std::endl;
  
  // Instantiate MiniDst reader
  MpdMiniDstReader *miniDstReader = new MpdMiniDstReader(inFileList);
  // Reader initialization
  miniDstReader->Init();
    // Specify branches to read
  miniDstReader->SetStatus("*", 0);
  miniDstReader->SetStatus("Event", 1);
  miniDstReader->SetStatus("Track", 1);
  miniDstReader->SetStatus("McTrack", 1);
  miniDstReader->SetStatus("BTofPidTraits", 1);

  if ( !miniDstReader->chain() ) {
    std::cout << "No chain has been found. Terminating..." << std::endl;
    return;
  }

  Long64_t events2read = miniDstReader->chain()->GetEntries();
  cout << "Number of events in chain: " << events2read << endl;

  //
  // MpdFemtoMaker setup
  //

  // MpdFemtoMaker
  MpdFemtoMaker* hbtMaker = new MpdFemtoMaker( "MpdFemto", "title" );

  // Create manager
  MpdFemtoManager* hbtManager = hbtMaker->hbtManager();

  // Setup MpdFemtoReader
  MpdFemtoMiniDstReader *miniHbtReader = new MpdFemtoMiniDstReader( miniDstReader );

  // Add reader to the manager
  hbtManager->setEventReader( miniHbtReader );

  // Create and set analysis
  MpdFemtoAnalysis *hbtAnalysis = new MpdFemtoAnalysis();
  hbtAnalysis->setVerboseMode( false );

  // Create and set event cut
  MpdFemtoBasicEventCut *eventCut = new MpdFemtoBasicEventCut();
  //eventCut->setVertZPos( -75., 75. );
  //eventCut->setEventMult( 0, 100 );
  hbtAnalysis->setEventCut( eventCut );

  // Create and set track cut for pion selection using TPC+TOF method
  MpdFemtoBasicTrackCut *trackCut = new MpdFemtoBasicTrackCut();
  trackCut->selectPrimary( true ); // analyse only primary tracks
  trackCut->setCharge( 1 );
  trackCut->setNHits( 15, 60 );
  trackCut->setEta( -1., 1. );
  trackCut->setPt( 0.15, 1.45 );
  trackCut->setP( 0.15, 1.45 );
  trackCut->setDCA( 0., 3. );
  trackCut->setDetectorSelection( 2 );
  trackCut->setTnTNSigmaPion( -2., 2. );
  trackCut->setMassSqr( -0.05, 0.08 );
  trackCut->setHbtPid( MpdFemtoBasicTrackCut::HbtPID::Pion );
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
  pairCut->setKt( 0.15, 0.66 );
  // Other antimerging/splitting and kinematic cuts

  // Create pair cut monitors
  MpdFemtoFxtPairCutMonitor *pairMoniPass = new MpdFemtoFxtPairCutMonitor("_moniPass_");
  MpdFemtoFxtPairCutMonitor *pairMoniFail = new MpdFemtoFxtPairCutMonitor("_moniFail_");
  pairCut->addCutMonitor(pairMoniPass, pairMoniFail);

  // Add pair cut to the analysis
  hbtAnalysis->setPairCut( pairCut );

  // Set how many events to mix
  hbtAnalysis->setNumEventsToMix( 5 );

  int kTBins = 5;
  double kTRange[2] = { 0.15, 0.65 };
  int qInvNbins = 50;
  double qInvRange[2] = { 0., 0.5 };
  MpdFemtoQinvCorrFctnKt *qInvCorrFctnKt =
    new MpdFemtoQinvCorrFctnKt( "qinv",
				qInvNbins, qInvRange[0], qInvRange[1],
				kTBins, kTRange[0], kTRange[1] );
  int qBPNbins = 40;
  double qBPRange[2] = { -0.2, 0.2 };
  MpdFemtoBPLCMS3DCorrFctnKt *qosl =
    new MpdFemtoBPLCMS3DCorrFctnKt( "bpLCMS", qBPNbins, qBPRange[0], qBPRange[1],
				    kTBins, kTRange[0], kTRange[1] );

  // Add correlation function to the analysis
  hbtAnalysis->addCorrFctn( qInvCorrFctnKt );
  hbtAnalysis->addCorrFctn( qosl );			      

  // Add analysis to the manager
  hbtManager->addAnalysis(hbtAnalysis);

  // Initialize MpdFemtoMaker
  int iret = hbtMaker->init();
  // Loop over events
  for (int iEvent = 0; iEvent < events2read; iEvent++) {
    
    hbtMaker->clear();
    std::cout << "hbtMiniDst -- Working on event # " << iEvent << std::endl;
    iret = hbtMaker->make();
    
    if (iret != 0) {
      std::cout << "[WARNING] MpdFemtoMaker returned status: " << iret
		<< ". Processing will be finished" << std::endl;
      break;
    }
  }

  // Maker finish flag
  iret = hbtMaker->Finish();

  // Write histograms to the output file
  TFile *oFile = new TFile( oFileName, "recreate" );
  trackMoniPass->writeOutHistos();
  trackMoniFail->writeOutHistos();
  pairMoniPass->writeOutHistos();
  pairMoniFail->writeOutHistos();
  qInvCorrFctnKt->writeOutHistos();
  qosl->writeOutHistos();
  oFile->Close();
}
