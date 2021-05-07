/*
 * MpdDstCompressTask2.h
 *
 *  Created on: 31 mar 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_NICA_HELPERS_MPDDSTCOMPRESSTASK2_H_
#define MPDROOT_NICA_MPD_NICA_HELPERS_MPDDSTCOMPRESSTASK2_H_
#include "FairTask.h"
#include "MpdEvent.h"
#include "FairEventHeader.h"
#include "FairLogger.h"
#include "MpdPid.h"
#include <TClonesArray.h>
class MpdDstCompressTask2 {
    Bool_t fUseMC;
    Bool_t fUseTpcKalmans;
    Bool_t fUseTpcHits;
    Bool_t fUseHeader;
    Bool_t fMCCompression;
    Bool_t fFixMom;
    MpdEvent *fMpdEvent;
    TClonesArray *fMCTracks;
    TClonesArray *fTpcKalmans;
    TClonesArray *fTpcHits;
    FairEventHeader *fEventHeader;
    Int_t fMCMapSize;
    Int_t *fMCIndexMap; //[fMCMapSize]
    MpdPid *fPID;
    virtual InitStatus CheckBranches();
    void Rewrite();
    void Fix();
    void FixPID(MpdTrack *track);
    void FixMom(MpdTrack* track, Int_t index, TVector3 &vertex);
public:
    MpdDstCompressTask2();
    void RegisterMC(Bool_t compress=kFALSE){fUseMC=kTRUE;fMCCompression=compress;};
    void RegisterTpcKalmans(){fUseTpcKalmans = kTRUE;};
    void RegisterTpcHits(){fUseTpcHits = kTRUE;};
    void RegisterEventHeader(){fUseHeader = kTRUE;};
    virtual InitStatus Init();
    virtual void  Exec(Option_t *option);
    virtual ~MpdDstCompressTask2();
};

#endif /* MPDROOT_NICA_MPD_NICA_HELPERS_MPDDSTCOMPRESSTASK2_H_ */
