//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      TpcClusterFinderAZ reads in TPC digits and reconstructs clusters and hits
//
//
// Environment:
//      Software developed for the MPD Detector at NICA.
//
// Author List:
//      Alexandr Zinchenko LHEP, JINR, Dubna - 14-July-2015
//
//-----------------------------------------------------------

#ifndef MPDTPCCLUSTERFINDERAZ_HH
#define MPDTPCCLUSTERFINDERAZ_HH

// Base Class Headers ----------------
#include "FairTask.h"
#include <TVector3.h>
#include <set>
#include <map>

// Collaborating Class Headers -------
#include "TpcSectorGeoAZ.h"

// Collaborating Class Declarations --
class MpdTpc2dCluster;
// class TpcGas;
class TClonesArray;

class MpdTpcClusterFinderAZ : public FairTask {
public:
   // Constructors/Destructors ---------
   MpdTpcClusterFinderAZ(BaseTpcSectorGeo &secGeo);
   ~MpdTpcClusterFinderAZ();

   // Operators

   // Accessors -----------------------

   // Modifiers -----------------------
   void SetPersistence(Bool_t opt = kTRUE) { fPersistence = opt; }

   // Operations ----------------------

   virtual InitStatus Init();
   void               FinishTask();

   virtual void Exec(Option_t *opt);
   // virtual void Clear(Option_t* opt);

private:
   // Private Data Members ------------
   static const Int_t fgkNsec2 = 24;                   // number of readout sectors (12 * 2)
   static const Int_t fgkNpads = 128, fgkNtimes = 512; // max number of pads and time bins

   TClonesArray *fDigiArray;
   TClonesArray *fClusArray;
   TClonesArray *fHitArray;
   // TClonesArray** fPrimArray;
   std::set<Int_t> *fDigiSet[fgkNsec2];
   Double_t         fCharges[fgkNpads][fgkNtimes];
   Int_t            fFlags[fgkNpads][fgkNtimes];
   Int_t            fDigis[fgkNpads][fgkNtimes];
   // TpcGas* fGas;
   Bool_t          fPersistence;
   TpcSectorGeoAZ *fSecGeo;

   // Private Methods -----------------
   void ProcessPadrow(Int_t isec, Int_t irow);                     // process one padrow of a sector
   void NextPixel(MpdTpc2dCluster *clus, Int_t ipad, Int_t itime); // add next pixel to the cluster
   void FindHits();                                                // find hits
   void FindHitsLocMax();                                          // find hits (local maxima)
   void PeakAndValley(const MpdTpc2dCluster *clus, std::multimap<Double_t, Int_t> &localMax); // peak-and-valley
   void CorrectReco(TVector3 &p3loc, TVector3 &p3err, Int_t nPads,
                    Double_t &adc); // correct reco charge, coordinates and errors

   ClassDef(MpdTpcClusterFinderAZ, 1);
};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
