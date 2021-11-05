/**
 * \brief Example of how to read a file (list of files) using McDst classes
 *
 * analyseMcDst.C is an example of reading McDst format.
 * One can use either uDst file or a list of mcDst files (inFile.lis or
 * inFile.list) as an input, and preform physics analysis
 */


// inFile - is a name of name.uDst.root file or a name
//          of a name.lis(t) files that contains a list of
//          name1.uDst.root files

#include <Rtypes.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

using namespace std;

R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "macro/mpd/mpdloadlibs.C"

//_________________
void analyseMcDst(const Char_t *inFile = "",
		  const Char_t *oFileName = "oProcTest.root") {

  McDstReader* myReader = new McDstReader(inFile);
  myReader->Init();

  // This is a way if you want to spead up IO
  cout << "Explicit read status for some branches" << endl;
  myReader->setStatus("*",0);
  myReader->setStatus("Event",1);
  myReader->setStatus("Particle",1);
  cout << "Status has been set" << endl;

  if( !myReader->chain() ) {
    cout << "No chain has been found." << endl;
  }
  Long64_t eventsInTree = myReader->tree()->GetEntries();
  cout << "eventsInTree: "  << eventsInTree << endl;
  Long64_t events2read = myReader->chain()->GetEntries();

  cout << "Number of events to read: " << events2read << endl;

  TFile *oFile = new TFile(oFileName, "RECREATE");

  // Histogramming
  // Event
  TH1D *hNch = new TH1D("hNCh","Number of charged particles;Nch;Entries",
                        300, -0.5, 2000.);

  // Track
  TH1D *hPz = new TH1D("hPz","p_{z} of particle;p_{z} (GeV/c);Entries",
                            402, -201., 201.);
  TH2D *hPtVsEta = new TH2D("hPtVsEta",
			    "p_{T} vs. #eta of primary track;#eta;p_{T} (GeV/c)",
                            220, -1.1, 1.1, 80, 0.05, 2.05);
  TH1D *hPionMom = new TH1D("hPionMom","Momentum of #pi;p (GeV/c);Entries",
                            100, 0., 2.);

  Int_t eventCounter = 0;
  Int_t hundredIter = 0;  

  // Loop over events
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {

    eventCounter++;
    if( eventCounter >= 100000 ) {
      eventCounter = 0;
      hundredIter++;
      cout << "Working on event #[" << (hundredIter * 100000)
        	      << "/" << events2read << "]" << endl;
    }

    Bool_t readEvent = myReader->loadEntry(iEvent);
    if( !readEvent ) {
      cout << "Something went wrong, Master! Nothing to analyze..."
								<< endl;
      break;
    }

    // Retrieve femtoDst
    MpdMcDst *dst = myReader->mcDst();

    // Retrieve event information
    McEvent *event = dst->event();
    if( !event ) {
      cout << "Something went wrong, Master! Event is hiding from me..."
		<< endl;
      break;
    }

    // Track analysis
    Int_t nTracks = dst->numberOfParticles();
    Int_t NumOfCharged = 0;
   
    cout << "Event# " << iEvent << " nTracks# " << nTracks << endl;
    
    // Track loop
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {        

      // Retrieve i-th femto track
      McParticle *particle = dst->particle(iTrk);

      if (!particle) continue;

      hPz->Fill( particle->pz() );
     
      if ( particle->charge() ) {
        NumOfCharged++;

        hPtVsEta->Fill( particle->eta(), particle->pt() );

        if ( particle->pdg() == 211 ) {
          hPionMom->Fill( particle->ptot() );
        }
      } // if ( particle->charge() )
    } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

    cout << NumOfCharged << endl;
    hNch->Fill( NumOfCharged );

  } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

  oFile->Write();
  oFile->Close();

  myReader->Finish();
}
