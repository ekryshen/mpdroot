
#ifndef MPD_VERTEXZFINDER_H
#define MPD_VERTEXZFINDER_H

#include "FairTask.h"

#include "TClonesArray.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TpcSectorGeoAZ.h"

class MpdVertexZfinder : public FairTask {
public:
   /** Constructor **/
   MpdVertexZfinder(BaseTpcGeo& fSecGeo, const char *name = "MpdVertexZfinder", Int_t iVerbose = 1);

   /** Destructor **/
   virtual ~MpdVertexZfinder();

   /// * FairTask methods

   /** Intialisation at begin of run. To be implemented in the derived class.
    *@return  Success   If not kSUCCESS, task will be set inactive.
    **/
   InitStatus Init();

   /** Reinitialisation.
    *@return  Success   If not kSUCCESS, task will be set inactive.
    **/
   InitStatus ReInit();

   /** Intialise parameter containers.
    **/
   void SetParContainers();

   void Exec(Option_t *option);

   /** Action after each event. **/
   void     Finish();
   void     Reset();
   void     SetHits(const TClonesArray *hits, const TH1F *hLays); // set hits container and histo
   Double_t FindZ(const Int_t *layPointers, Int_t &flag);         // evaluate vertex Z-position

private:
   TpcSectorGeoAZ     *secGeo;  // tpc sector geometry
   const TClonesArray *fKHits;  // array of Kalman hits
   const TH1F         *fhLays;  // histogram of layer occupancy
   TH1F               *fhZ;     // histogram of Z-positions
   TH2F               *fhZDip;  // histogram Z - Dip angle
   TF1                *fUnc;    // fitting function
   TH2F               *fhPhLay; // histogram of phase space occupancy

private:
   // Some constants
   // static const Double_t fgkChi2Cut; // max accepted Chi2 of hit for track

   ClassDef(MpdVertexZfinder, 1);
};
#endif
