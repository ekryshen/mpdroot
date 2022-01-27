/*
 * pt_spec.C
 *
 *  Created on: 27 sty 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 *
 *
 *      Example of FairTask usage for simple analysis
 */

#ifndef __CLING__
#include <FairRunAna.h>
#include <vector>
#include <RtypesCore.h>
#include <TString.h>
#include "MpdMiniDstSource.h"
#include <FairTask.h>
#include <TH1D.h>
#include "../../../physics/common/v0/MpdV0CandidateCut.h"
#include "../../../physics/common/v0/MpdV0FinderBasic.h"
#include "../../../physics/common/minidst/MpdMiniDstSource.h"
#include "../../../physics/common/v0/MpdV0DaughterCutBasic.h"
#include "../../../physics/common/v0/MpdV0FinderHelix.h"
#include "../../../physics/common/v0/MpdV0Namespace.h"
#include "../../../physics/nicafemto/nica_helpers/MpdPIDOnTheFly.h"
#include "../../../external/nicafemto/features/NicaStd.h"

#endif

class SimpleSpectra : public FairTask {
   TH1D *        fPt;
   TClonesArray *fMiniDstTracks;

protected:
   /** this method is called when you call FairRunAna::Init in principle here you are "connecting" to the data,
    * initialize histograms
    *
    */
   virtual InitStatus Init()
   {
      fMiniDstTracks = (TClonesArray *)FairRootManager::Instance()->GetObject("Track");
      if (fMiniDstTracks == nullptr) return kFATAL;
      fPt = new TH1D("pt", "pt", 100, 0, 10);
      return kSUCCESS;
   };
   /**
    * method called at the end of analysis, here you should write the output of your analysis
    */
   virtual void Finish() { fPt->Write(); };

public:
   /**
    * Default constructor, remember to set all pointers to null if you are planning to use proof!
    */
   SimpleSpectra() : fPt(nullptr), fMiniDstTracks(nullptr){};
   /**
    * this method is called for each entry in processed tree, here you are doing the analysis
    */
   virtual void Exec(Option_t *option)
   {
      for (int i = 0; i < fMiniDstTracks->GetEntriesFast(); i++) {
         MpdMiniTrack *track = (MpdMiniTrack *)fMiniDstTracks->UncheckedAt(i);
         fPt->Fill(track->gPt());
      }
   }
   virtual ~SimpleSpectra(){};
};
void pt_spec(TString inFile  = "/media/daniel/hdd/minidst/urqmd-AuAu-07.7GeV-mb-eos0-500-960.reco.MiniDst.root",
             TString outFile = "spectra.root")
{
   FairRunAna *run = new FairRunAna();
   run->SetEventHeaderPersistence(kFALSE); // don't write event header

   MpdMiniDstSource *source = new MpdMiniDstSource(inFile);
   run->SetSource(source);
   run->SetOutputFile(outFile);
   run->AddTask(new SimpleSpectra());
   run->Init();
   run->Run();
}
