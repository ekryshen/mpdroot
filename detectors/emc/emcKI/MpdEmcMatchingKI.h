//--------------------------------------------------------------------
//
// Description:
//      MPD TPC-EMC Matching
//
//
// Author List: D.Peresunko
//
//--------------------------------------------------------------------

#ifndef MPDEMCMATCHINGKI_H
#define MPDEMCMATCHINGKI_H 1

#include "FairTask.h"

#include <vector>

#include "MpdEmcTrackExtrap.h"
class MpdEmcGeoParams;
class MpdTpcKalmanTrack;
class TClonesArray;

class MpdEmcMatchingKI : public FairTask
{
 public:
  MpdEmcMatchingKI();

  ~MpdEmcMatchingKI();

  virtual InitStatus Init();

  virtual void Exec(Option_t* opt);
  void virtual Finish();

 protected:
  void ExtrapolateTracks();   // Fill array with points of track extapolations to ECAL surface
  void MakeMatchToClusters(); // find best match to clusters
  void MakeMatchToTracks();   // Find best match to tracks
  void CorrectClustersVtx();  // correct clulster position for Zcoordinate of vertex

 private:
  std::vector<MpdEmcTrackExtrap> fTrackPoints; //! Coordinates and link to track extrapolations
  TObjArray* fClusterArray;                    //! Input/output list of clusters
  TClonesArray* fTpcTracks;                    //! Input/output clusters array
  TClonesArray* fvtx;                          //! list of reconstructed vertexes
  double fzMax;                                // Maximal size of ECAL in z direction
  double frMin;                                // inner radius of ECAL

  ClassDef(MpdEmcMatchingKI, 1);
};

#endif
