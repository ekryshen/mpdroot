#ifndef MPD_KFPRIMARYVERTEXFINDER_H
#define MPD_KFPRIMARYVERTEXFINDER_H

#include "MpdKalmanTrack.h"

#include "FairTask.h"

#include "TMatrixD.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TH1D.h"

class MpdKfPrimaryVertexFinder : public FairTask {
public:
   /** Constructor **/
   MpdKfPrimaryVertexFinder(const char *name = "MpdKfPrimaryVertexFinder", Int_t iVerbose = 1);

   /** Destructor **/
   virtual ~MpdKfPrimaryVertexFinder();

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
   void Finish();
   void Reset();

   // Set containers - to use the task as an "external" one
   void SetVertices(TClonesArray *vertCont) { fVertexCont = vertCont; } // set vertex container
   void SetTracks(TClonesArray *tracks) { fTracks = tracks; }           // set track container
   void SetConstrFlag(Int_t flag = 1) { fConstrFlag = flag; }           // set vertex constrain flag
   void SetSmoothSame(Int_t flag = 1) { fSmoothSame = flag; }           // smooth original tracks
   void Smooth();     // smooth (update track momenta and track lengths)
   void Chi2Vertex(); // compute chi2-distance of tracks to primary vertex

private:
   void EvalVertex(); // evaluate primary vertex position
   void FindVertex(); // find primary vertex
   void ComputeAandB(TMatrixD &xk0, const MpdKalmanTrack *track, const MpdKalmanTrack &trackM, TMatrixD &a, TMatrixD &b,
                     TMatrixD &ck0);                                 // compute matrices of derivatives
   void Proxim(const MpdKalmanTrack &track0, MpdKalmanTrack &track); // adjust track params
   void FillVertex();                                                // fill vertex info
   void Write();
   void Writedir2current(TObject *obj);

   Int_t       fConstrFlag; // flag for storing tracks with vertex constrain ( if != 0)
   Int_t       fSmoothSame; // flag for smoothing original tracks (if fConstrFlag == 0)
   Int_t       fNPass;      //
   Double_t    fXYZ[3];     // vertex position
   TMatrixD    fCovar;      // vertex position covariance matrix
   TF1        *fUnc;        // fitting function
   TH1D       *fHist[3];    // X-,Y-,Z-histograms
   TDirectory *fHistoDir;   //

   TClonesArray *fTracks;     // TPC tracks
   TClonesArray *fMCTracks;   // MC tracks
   TClonesArray *fVertexCont; // vertices
   TClonesArray *fVertTracks; // tracks with vertex constrain
                              // TClonesArray *fSTSTrackMatch;

   // FairVertex *fPrimVtx;
   // FairStsKFTrackFitter fStsFitter;

private:
   // Some constants
   static const Double_t fgkChi2Cut; // max accepted Chi2 of track to be added to the vertex

   ClassDef(MpdKfPrimaryVertexFinder, 2);
};

#endif
