// skelet

#ifndef ROOT_MpdRoInvMassTask
#define ROOT_MpdRoInvMassTask
#ifndef ROOT_FairTask
#include "FairTask.h"
#endif

#include <TDatabasePDG.h>
#include <TH1F.h>

#include "MpdEvent.h"

class MpdRoInvMassTask : public FairTask {

private:
    UInt_t fEventCounter; //! event counter
    TDatabasePDG*  fPDG;  //!  PDG database

    MpdEvent *fDstEvent; // dst event
    TClonesArray *fMCTracks; // array of MC tracks
    
    TH1F *fRoInvMass;
    TH1F *fRoInvMassMC;

public:

  MpdRoInvMassTask();
  MpdRoInvMassTask(const char *name, const char *title="MPD Analysis");
  virtual ~MpdRoInvMassTask();	// Destructor

    virtual void Exec(Option_t * option);
    virtual InitStatus Init(); // Init before Exec
    virtual void Finish(); // Finish after Exec

    void  Reset(); //
    void  Register(); // Register what do you want to write to file
    void  SetOption(Option_t *option=" ") {fOption = option;  fOption.ToLower();}

  ClassDef(MpdRoInvMassTask,0)
};

#endif
