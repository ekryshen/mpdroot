/*
 * MpdDstWriteTask.h
 *
 *  Created on: 23 lut 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDDST_MPDDSTCOMPRESSTASK_H_
#define MPDDST_MPDDSTCOMPRESSTASK_H_

#include "FairTask.h"
#include "MpdEvent.h"
#include "FairEventHeader.h"
#include "FairLogger.h"
#include <TClonesArray.h>

/**
 * class to write only miniaml output to new root file
 * (MpdEvent + event header)
 */
class MpdDstCompressTask : public FairTask {
protected:
   Bool_t             fUseMC;
   Bool_t             fUseFreezouts;
   Bool_t             fUseTpcKalmans;
   Bool_t             fUseTpcHits;
   Bool_t             fUseHeader;
   Bool_t             fMCCompression;
   MpdEvent          *fMpdEvent;
   TClonesArray      *fFreezouts;
   TClonesArray      *fMCTracks;
   TClonesArray      *fTpcKalmans;
   TClonesArray      *fTpcHits;
   Int_t              fMCMapSize;
   Int_t             *fMCIndexMap; //[fMCMapSize]
   virtual InitStatus CheckBranches();

public:
   MpdDstCompressTask();
   MpdDstCompressTask(const char *name, Int_t Verbose = 1);
   void RegisterMC(Bool_t compress = kFALSE)
   {
      fUseMC         = kTRUE;
      fMCCompression = compress;
   };
   void               RegisterMCFreezouts() { fUseFreezouts = kTRUE; };
   void               RegisterTpcKalmans() { fUseTpcKalmans = kTRUE; };
   void               RegisterTpcHits() { fUseTpcHits = kTRUE; };
   void               RegisterEventHeader() { fUseHeader = kTRUE; };
   virtual InitStatus Init();
   virtual void       Exec(Option_t *option);
   virtual ~MpdDstCompressTask();
   ClassDef(MpdDstCompressTask, 1);
};

#endif /* MPDDST_MPDDSTCOMPRESSTASK_H_ */
