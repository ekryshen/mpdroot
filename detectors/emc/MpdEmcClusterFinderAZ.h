//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      EmcClusterFinderAZ reads in EMC digits and reconstructs clusters
//
//
// Environment:
//      Software developed for the MPD Detector at NICA.
//
// Author List:
//      Alexandr Zinchenko LHEP, JINR, Dubna - 18-May-2016
//      Alexandr Zinchenko LHEP, JINR, Dubna - 24-June-2018 - adapted for projective geometry 
//
//-----------------------------------------------------------

#ifndef MPDEMCCLUSTERFINDERAZ_HH
#define MPDEMCCLUSTERFINDERAZ_HH

// Base Class Headers ----------------
#include "FairTask.h"
#include <TVector3.h>
#include <set>
#include <map>

// Collaborating Class Headers -------

// Collaborating Class Declarations --
class MpdTpc2dCluster;
class MpdEmcGeoParams;
class TClonesArray;

class MpdEmcClusterFinderAZ : public FairTask {
public:

  // Constructors/Destructors ---------
  MpdEmcClusterFinderAZ();
  ~MpdEmcClusterFinderAZ();

  // Operators
  
  // Accessors -----------------------

  // Modifiers -----------------------
  void SetPersistence(Bool_t opt = kTRUE) { fPersistence = opt; }
  void SetThreshold(Double_t thresh) { fThresh = thresh; }
  void SetTimeWindow(Int_t iwind = 1) { fTimeWindow = iwind; }

  // Operations ----------------------
  
  virtual InitStatus Init();
  void FinishTask();

  virtual void Exec(Option_t* opt);
  //virtual void Clear(Option_t* opt);

private:

  TClonesArray* fDigiArray;
  TClonesArray* fClusArray;
  TClonesArray* fHitArray;
  TClonesArray* fMcTrArray;
  //TClonesArray** fPrimArray;
  /*
  std::set<Int_t>* fDigiSet[fgkNsec2];
  Double_t fCharges[fgkNpads][fgkNtimes];
  Int_t fFlags[fgkNpads][fgkNtimes];
  Int_t fDigis[fgkNpads][fgkNtimes];
  */
  std::vector<std::set<Int_t> > fDigiSet;
  std::vector<std::vector<Double_t> > fCharges;
  std::vector<std::vector<Int_t> > fFlags, fDigis;
  std::vector<Int_t> fSecRows0;
  std::set<Int_t> fSecRows1;
  Bool_t fPersistence;
  MpdEmcGeoParams* fEmcGeo;
  Float_t fThresh;
  Int_t fTimeWindow;

  // Private Methods -----------------
  void FillEmcInfo(); // fill EMC info
  void ProcessSector(Int_t isec); // process one sector
  void NextPixel(MpdTpc2dCluster* clus, Int_t iphi, Int_t iz); // add next pixel to the cluster
  void FindHits(); // find hits
  void GetPhiTheta(Double_t &phi, Double_t &theta); // convert COG in bin units to angles
  void FindHitsLocMax(); // find hits (local maxima)
  void PeakAndValley(const MpdTpc2dCluster* clus, std::multimap<Double_t,Int_t> &localMax, Int_t ishift); // peak-and-valley
  void RedoId(std::map<Int_t,Float_t>& contrib); // recompute ID contributions (for tracks born outside EMC) 

  ClassDef(MpdEmcClusterFinderAZ,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
