#ifndef MPD_TPCDEDXTASK_H
#define MPD_TPCDEDXTASK_H

#include "FairTask.h"

#include "TClonesArray.h"

class MpdTpcKalmanTrack;

class MpdTpcDedxTask : public FairTask {
public:
   /** Constructor **/
   MpdTpcDedxTask(const char *name = "MpdTpcDedxTask", Int_t iVerbose = 1);

   /** Destructor **/
   virtual ~MpdTpcDedxTask();

   /// * FairTask methods

   /** Intialisation at begin of run. To be implemented in the derived class.
    *@value  Success   If not kSUCCESS, task will be set inactive.
    **/
   InitStatus Init();

   /** Intialise parameter containers.
    **/
   void SetParContainers();

   void Exec(Option_t *option);

   /** Action after each event. **/
   void Finish();
   void Reset();

private:
   void DriftCorrection(MpdTpcKalmanTrack *track); // hit correction for drift length
   void Write();
   void Writedir2current(TObject *obj);

   TDirectory *fHistoDir; // pointer to histogram directory

   TClonesArray *fTracks;   // TPC tracks
   TClonesArray *fMCTracks; // MC tracks
   TClonesArray *fHits0;    // TPC hits (from hit producer)
   TClonesArray *fHits;     // TPC hits (rec. points)

private:
   ClassDef(MpdTpcDedxTask, 0);
};

#endif
