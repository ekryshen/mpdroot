// -------------------------------------------------------------------------
// -----                 MpdEmcClusterCreator header file                 -----
// -------------------------------------------------------------------------

#ifndef MPDEMCCLUSTERCREATION_H
#define MPDEMCCLUSTERCREATION_H 1

#include <iostream>
#include "FairTask.h"
#include "MpdEmcGeoParams.h"

class TClonesArray;

class MpdEmcClusterCreation : public FairTask {
public:
   /** Default constructor **/
   MpdEmcClusterCreation();

   /** Destructor **/
   ~MpdEmcClusterCreation();

   /** Virtual method Init **/
   virtual InitStatus Init();

   /** Virtual method Exec **/
   virtual void Exec(Option_t *opt);
   void virtual Finish();

   // Search relative module of the central hit inside frame with sides rowFrame and modFrame
   void SearchFrameHits(Int_t row, Int_t mod, vector<Int_t> &relHits);

   // Search relative module of the central hit with row/mod position
   void SearchRelativeHits(Int_t row, Int_t mod, vector<Int_t> &relHits);

   // Set threshold for each hit
   void SetEnergyThreshold(Float_t fEnMin) { fEnergyThreshold = fEnMin; };

   // Set max radius for each cluster
   void SetMaxClusterRadius(Float_t fRad) { fMaxClusterRadius = fRad; };

   // Set frame to search cluster candidates
   void SetClusterFrame(Float_t row, Float_t mod)
   {
      rowFrame = row;
      modFrame = mod;
   };

   Float_t GetEnergyThreshold() { return fEnergyThreshold; };

   Float_t GetMaxClusterRadius() { return fMaxClusterRadius; };

private:
   /** Input array of MpdEmcHit **/
   TClonesArray *fHitArray;
   /** Input array of MCTracks **/
   TClonesArray *fMcTrackArray;
   /** Output array of MpdEmcCluster **/
   TClonesArray *fClusterArray;

   Float_t fEnergyThreshold; // Energy threshold for each module

   // First method variable
   Float_t fMaxClusterRadius; // Maximal radius of cluster

   // Second method variables
   UInt_t rowFrame; // row window to search cluster candidates
   UInt_t modFrame; // tower (z) window to search cluster candidates

   MpdEmcGeoParams *fGeoPar;

   ClassDef(MpdEmcClusterCreation, 2);
};

#endif
