// C++ headers
#include <vector>
// ROOT headers
#include "Rtypes.h"
#include "TChain.h"
#include "TFile.h"
// MpdRoot part
R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "macro/mpd/mpdloadlibs.C"

//
// inFineName can take either filename.MiniDst.root file or
// a list of such files with .list or .lis extention (for example,
// somelist.list or somelist.lis)
//

//_________________
void miniDstProcExample(const Char_t* inFileName) {
  
  // Instantiate reader
  MpdMiniDstReader* miniDstReader = new MpdMiniDstReader(inFileName);
  
  // Reader initialization
  miniDstReader->Init();
  
  // One can specify branches to read
  miniDstReader->SetStatus("*", 0);               // Turn off all branches
  miniDstReader->SetStatus("Event*", 1);          // Turn on specific branch
  miniDstReader->SetStatus("Track*", 1);
  miniDstReader->SetStatus("BTofHit*", 1);
  miniDstReader->SetStatus("BTofPidTraits*", 1);
  miniDstReader->SetStatus("BECalCluster*", 1);
  miniDstReader->SetStatus("FHCalHit*", 1);
  miniDstReader->SetStatus("TrackCovMatrix*", 0); // Turn off specific branch
  miniDstReader->SetStatus("McEvent*", 1);
  miniDstReader->SetStatus("McTrack*", 1);
  
  // Retrieve events in tree
  Long64_t events2read = miniDstReader->chain()->GetEntries();

  // Loop over events
  for (Long64_t i = 0; i < events2read; i++) {

    // Read next event
    Bool_t  isOk = miniDstReader->readMiniEvent( i );
    
    // Retrieve current miniDst (from the given .MiniDst.root file)
    MpdMiniDst *dst = miniDstReader->miniDst();

    //
    // Event information
    //
    
    // Get MiniEvent
    MpdMiniEvent *event = dst->event();

    // Get primary vertex z-position
    Float_t z = event->primaryVertex().Z();

    //
    // Track information
    //

    // Retrieve number of reconstructed tracks
    Int_t nGTracks = dst->numberOfTracks();
    
    // Track loop
    for (Int_t j = 0; j < nGTracks; j++) {
      
      // Retrieve j-th mini track
      MpdMiniTrack *miniTrack = dst->track(j);
      
      // Global track full momentum
      Float_t ptot = miniTrack->gMom().Mag();

      //
      // TOF-matched tracks
      //

      // Check if track matched TOF
      if ( miniTrack->isBTofTrack() ) {
	
	// Retrieve TOF information for the matched track
	MpdMiniBTofPidTraits *tofTrait = dst->btofPidTraits( miniTrack->bTofPidTraitsIndex() );
	
	// Retrive beta
	Float_t beta = tofTrait->beta();
      }
    }


    //
    // Hit information
    //
    
    // Loop over barrel TOF hits
    for (Int_t j = 0; j < dst->numberOfBTofHits(); j++ ) {
      
      // Retrieve j-th hit information
      MpdMiniBTofHit *btofHit = dst->btofHit( j );
      
      // x position of hit
      Float_t xPosition = btofHit->btofHitPosX();
    }
    
    // Loop over FHCal hits
    for (Int_t j = 0; j < dst->numberOfFHCalHits(); j++) {

      // Retrieve j-th hit information
      MpdMiniFHCalHit *hit = dst->fhcalHit( j );

      // Module number that corresponds to the hit
      Int_t modId = hit->module();
    }

    //
    // BECal cluster information
    //

    // Loop over clusters
    for (Int_t j = 0; j < dst->numberOfBECalClusters(); j++) {

      // Retrieve j-th cluster
      MpdMiniBECalCluster *cluster = dst->becalCluster( j );

      // Retrieve cluster energy
      Float_t energy = cluster->energy();

      // Retrieve mcTrack IDs that produced cluster
      std::vector<UShort_t> mcTrackIds = cluster->mcTrackIds();
    }

    //
    // Monte Carlo information
    //

    // Retrieve MC event
    MpdMiniMcEvent *mcEvent = dst->mcEvent();
    // Impact parameter
    Float_t b = mcEvent->primaryVertexZ();

    // MC track loop
    for (Int_t j = 0; j < dst->numberOfMcTracks(); j++) {
      
      // Retrieve j-th MC track
      MpdMiniMcTrack *mcTrack = dst->mcTrack( j );

      // Print MC track information
      mcTrack->Print();
    }
    
  }
  
  // Finalize miniDst reader
  miniDstReader->Finish();
}
