//--------------------------------------------------------------------
//
// Description:
//      MPD TPC-EMC Matching 
//
//
// Author List:
//      Alexander Zinchenko LHEP, JINR, Dubna - 1-June-2016
//      Alexander Zinchenko LHEP, JINR, Dubna - 8-June-2018 - adapted for projective geometry
//
//--------------------------------------------------------------------

#ifndef MPDEMCMATCHING_H
#define MPDEMCMATCHING_H 1

//#include "MpdEmcDigit.h"

#include "FairTask.h"

#include <TVector3.h>
#include <iostream>
#include <map>
#include <set>

class MpdEmcGeoParams;
class MpdTpcKalmanTrack;
class TClonesArray;

class MpdEmcMatching : public FairTask 
{
 public:

  /** Default constructor **/
  MpdEmcMatching();
  
  
  /** Destructor **/
  ~MpdEmcMatching();
  
  
  /** Virtual method Init **/
  virtual InitStatus Init();
  
  
  /** Virtual method Exec **/
  virtual void Exec(Option_t* opt);
  void virtual Finish();
  
 private:

  /** Input array of MpdEmcPoints **/
  TClonesArray* fPointArray;
  /** Input array of MCTracks **/
  TClonesArray* fMcTrackArray;
  /** Input array of EmcRecPoints **/
  TClonesArray* fRecPointArray;
  /** Input array of TpcKalmanTracks **/
  TClonesArray* fTpcTracks;
  /** Output array of EmcMatches **/
  TClonesArray* fMatchArray;
  
  /** Output array of MpdEmcHit **/
  TClonesArray* fDigiArray;
  
  void DoMatching(Int_t itrack);
  void GetTowerCoords(TVector3 &pos, Int_t io, Double_t &phiT, Double_t &theT);
  
  MpdEmcGeoParams* fGeoPar;
  std::vector<std::multimap<Double_t,Int_t> > fRecPoints;
  std::set<Int_t> fSecRows0;
  
  ClassDef(MpdEmcMatching, 1);
  
};

#endif
