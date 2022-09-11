//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      TpcClusterFinderMlem reads in TPC digits and reconstructs clusters and hits
//
//
// Environment:
//      Software developed for the MPD Detector at NICA.
//
// Author List:
//      Alexandr Zinchenko LHEP, JINR, Dubna - 11-January-2016
//
//-----------------------------------------------------------

#ifndef MPDTPCCLUSTERFINDERMLEM_HH
#define MPDTPCCLUSTERFINDERMLEM_HH

// Base Class Headers ----------------
#include "FairTask.h"
#include <TMatrixD.h>
#include <TVector3.h>
#include <set>
#include <vector>
#include <map>

// Collaborating Class Headers -------
#include "TpcSectorGeoAZ.h"

// Collaborating Class Declarations --
class MpdTpc2dCluster;
class TClonesArray;
class TH2D;
// class TMatrixD;

class MpdTpcClusterFinderMlem : public FairTask {
public:
   // Constructors/Destructors ---------
   MpdTpcClusterFinderMlem(BaseTpcGeo& secGeo);
   ~MpdTpcClusterFinderMlem();

   // Operators

   // Accessors -----------------------

   // Modifiers -----------------------
   void SetPersistence(Bool_t opt = kTRUE) { fPersistence = opt; }

   // Operations ----------------------

   virtual InitStatus Init();
   void               FinishTask();

   virtual void Exec(Option_t *opt);
   // virtual void Clear(Option_t* opt);

   struct pixel {
      Double_t qq;
      Int_t    ix;
      Int_t    iy;
      Double_t vis;
      Double_t sum;
   };

private:
   // Private Data Members ------------
   static const Int_t fgkNsec2 = 24;                   // number of readout sectors (12 * 2)
   static const Int_t fgkNpads = 128, fgkNtimes = 512; // max number of pads and time bins
   // AZ static const Int_t fgkOvfw = 4095; // overflow value
   static const Int_t fgkOvfw = 1023; // overflow value - 10 bits

   TClonesArray *fDigiArray;
   TClonesArray *fClusArray;
   TClonesArray *fHitArray;
   // TClonesArray** fPrimArray;
   std::set<Int_t> *fDigiSet[fgkNsec2];
   Double_t         fCharges[fgkNpads][fgkNtimes];
   Int_t            fFlags[fgkNpads][fgkNtimes];
   Int_t            fDigis[fgkNpads][fgkNtimes];
   Bool_t           fPersistence;
   TpcSectorGeoAZ  *fSecGeo;
   Double_t         fVertexZ; // !!! vertex Z-position estimate (at present it is taken from true value) !!!

   // Private Methods -----------------
   void ProcessPadrow(Int_t isec, Int_t irow);                       // process one padrow of a sector
   void NextPixel(MpdTpc2dCluster *clus, Int_t ipad, Int_t itime);   // add next pixel to the cluster
   void FindHits();                                                  // find hits
   void Mlem(Int_t iclus, std::multimap<Double_t, Int_t> &localMax); // MLEM procedure
   // bool myfunc(const pixel i, const pixel j); // sorting function
   // static Bool_t SortPix(const pixel i, const pixel j); // sorting function
   void     GetResponse(const MpdTpc2dCluster *clus, TH2D *hXY, TH2D *hOvfw, Double_t &sigt, Double_t &sigp,
                        Double_t &correl); // get response parameters
   Double_t GetCij(Int_t irow, Double_t x0, Double_t y0, Double_t x1, Double_t y1, Double_t sigt, Double_t sigp,
                   Double_t correl, Double_t &vis); // compute pixel-to-bin couplings
   void     PeakAndValley(const MpdTpc2dCluster *clus, std::multimap<Double_t, Int_t> &localMax); // peak-and-valley
   void     PeakAndValley(const std::vector<pixel> &pixels, std::multimap<Double_t, Int_t> &localMax,
                          std::vector<std::vector<Double_t>> &charges,
                          std::vector<std::vector<Int_t>>    &flags); // peak-and-valley in pixel domain
   void     CreateHits(const std::vector<pixel> &pixels, std::multimap<Double_t, Int_t> &localMax,
                       std::vector<std::vector<Double_t>> &charges, std::vector<std::vector<Int_t>> &flags, Int_t iclus,
                       std::vector<std::multimap<Double_t, Int_t>> &pixInMax);    // create hits from pixels
   void CorrectReco(TVector3 &p3loc, TVector3 &p3err, Int_t nPads, Double_t adc); // correct reco coordinates and errors
   void CorrectRecoMlem(TVector3 &p3loc, TVector3 &p3errCor, MpdTpc2dCluster *clus, Double_t adc); // after MLEM
   void ChargeMlem(Int_t nHits0, std::vector<pixel> &pixels, std::vector<pixel> &bins,
                   std::vector<std::multimap<Double_t, Int_t>> &pixInMax, const TMatrixD &cij,
                   Double_t cijMin); // correct hit charges after MLEM

   ClassDef(MpdTpcClusterFinderMlem, 0);
};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------
