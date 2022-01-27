/*************************************************************************************
 *
 *         Class MpdZdcTstSim
 *
 *  Author:   Elena Litvinenko
 *  e-mail:   litvin@nf.jinr.ru
 *  Version:  8-Apr-2008
 *
 ************************************************************************************/

#ifndef MPDZDCTSTSIM_H
#define MPDZDCTSTSIM_H

#include "TClonesArray.h"
#include "TNtuple.h"
#include "TNtupleD.h"
#include "TH3F.h"

#include "FairTask.h"
#include "FairMCEventHeader.h"
#include "MpdZdcGeoPar.h"

class MpdZdcTstSim : public FairTask {

public:
   MpdZdcTstSim();
   MpdZdcTstSim(const char *name, const char *title, Int_t verbose);
   MpdZdcTstSim(const char *name, const char *title, Int_t verbose, Int_t flag);
   virtual ~MpdZdcTstSim();

   virtual InitStatus Init();
   virtual InitStatus ReInit();

   virtual void Exec(Option_t *option = "0");
   virtual void Finish();

   void CreateMyTree();

private:
   Int_t nevents;

   TClonesArray      *fMCTrackArray;    //!
   TClonesArray      *fMCZdcPointArray; //!
   TClonesArray      *fMCTpcPointArray; //!
   FairMCEventHeader *fMCEventHeader;   //!

   TNtuple  *fTree;        //
   TNtupleD *fTreeSummary; //
   TH3F     *fH3;          //

   //  MpdZdcGeoPar *fPar;           //  !

   ClassDef(MpdZdcTstSim, 1);
};
#endif
