/*
 * MpdV0Matcher.h
 *
 *  Created on: 27 sty 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_PHYSICS_COMMON_V0_MPDV0MATCHER_H_
#define MPDROOT_PHYSICS_COMMON_V0_MPDV0MATCHER_H_

#include <FairTask.h>
#include <Rtypes.h>
#include <RtypesCore.h>
#include <vector>

#include "TVector3.h"

#include "MpdSimpleLinks.h"
#include "MpdV0Namespace.h"

/**
 * class for matchin found v0 particles with MC tracks
 */
class MpdEvent;

class MpdV0Matcher : public FairTask {

public:
   enum class EMatchType { kMatchByMomentum, kMatchByFirstDaughter, kMatchBySecondDaughet, kMatchByBothDaughters };

protected:
   Bool_t                     fWrite;
   Bool_t                     fFirst;
   MpdCommonV0::EParticleType fPidHipo;
   EMatchType                 fMethod;
   Double_t                   fMomThreshold;
   MpdEvent *                 fMpdEvent;
   TClonesArray *             fMiniEvents;
   TClonesArray *             fMiniTracks;
   TClonesArray *             fMcTracks;
   TClonesArray *             fV0s;
   MpdSimpleLinks<Int_t> *    fLinks;
   std::vector<int>           fMcPdg;
   std::vector<int>           fMcV0Indexes;
   std::vector<int>           fMcV0sPid;
   std::vector<TVector3>      fMcV0Momenta;

   enum class EFormat { kDst, kMiniDst };
   EFormat            fFormat;
   Bool_t             IsGoodV0(Int_t pid) const;
   void               PrepareMomentumMatchingDst();
   void               PrepareMomentumMatchinMiniDst();
   void               MatchV0ByMomentum();
   void               MatchDstByDaugher(Int_t id);
   void               MatchMiniDstByDaughter(Int_t id);
   void               MatchDstByDaughers();
   void               MatchMiniDstByDaughters();
   void               MatchByMc();
   void               MatchByMomentum();
   virtual InitStatus Init();

public:
   MpdV0Matcher(MpdCommonV0::EParticleType pidHypo = MpdCommonV0::EParticleType::kPdgHypo,
                EMatchType                 type    = EMatchType::kMatchByMomentum);
   void SetPidHypo(MpdCommonV0::EParticleType pidHypo) { fPidHipo = pidHypo; };
   void SetMethod(EMatchType type) { fMethod = type; };
   void Write(Bool_t write) { fWrite = write; }
   /**
    * set maximum momentum difference (in percent) between found V0 and matched MC V0, if this value is bigger
    * than threshold v0 is marked as not matched
    */
   void         SetMomThreshold(Double_t threshold) { fMomThreshold = threshold; }
   virtual void Exec(Option_t *option);
   virtual ~MpdV0Matcher();
   ClassDef(MpdV0Matcher, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_MPDV0MATCHER_H_ */
