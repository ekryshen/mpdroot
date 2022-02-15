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
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TVector3.h>
#include <type_traits>
#include <vector>

#include "MpdMCTrack.h"
#include "MpdEvent.h"
#include "MpdTrack.h"
#include "MpdMiniMcTrack.h"
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
   Bool_t                 fWrite;
   Bool_t                 fFirst;
   MpdV0::EParticleType   fPidHipo;
   EMatchType             fMethod;
   Double_t               fMomThreshold;
   MpdEvent *             fMpdEvent;
   TClonesArray *         fMiniEvents;
   TClonesArray *         fMiniTracks;
   TClonesArray *         fMiniMcTracks;
   TClonesArray *         fMcTracks;
   TClonesArray *         fV0s;
   MpdSimpleLinks<Int_t> *fLinks;
   std::vector<int>       fMcPdg;
   std::vector<int>       fMcV0Indexes;
   std::vector<int>       fMcV0sPid;
   std::vector<TVector3>  fMcV0Momenta;

   enum class EFormat { kDst, kMiniDst };
   EFormat fFormat;
   Bool_t  IsGoodV0(Int_t pid) const;
   void    PrepareMomentumMatchingDst();
   void    PrepareMomentumMatchinMiniDst();
   /**
    * we have to formats, to not write everything twice let's use templates
    */
   template <class T>
   Int_t GetMotherIndex(T *val)
   {
      if constexpr (std::is_same<T, MpdMCTrack>::value) {
         return val->GetMotherId();
      } else {
         return val->getMotherId();
      }
   }
   template <class T>
   T *GetRecoTrack(Int_t index)
   {
      if constexpr (std::is_same<T, MpdTrack>::value) {
         return (T *)fMpdEvent->GetGlobalTracks()->UncheckedAt(index);
      } else {
         return (T *)fMiniTracks->UncheckedAt(index);
      }
   }
   template <class T>
   T *GetSimTrack(Int_t index)
   {
      if constexpr (std::is_same<T, MpdMCTrack>::value) {
         return (T *)fMcTracks->UncheckedAt(index);
      } else {
         return (T *)fMiniMcTracks->UncheckedAt(index);
      }
   }
   template <class T>
   inline Int_t GetPdgCode(const T *val)
   {
      if constexpr (std::is_same<T, MpdMCTrack>::value) {
         return val->GetPdgCode();
      } else {
         return val->pdgId();
      }
   }
   template <class T>
   inline Int_t GetMother(const T *val)
   {
      if constexpr (std::is_same<T, MpdMCTrack>::value) {
         return val->GetMotherId();
      } else {
         return val->getMotherId();
      }
   }
   template <class T1, class T2>
   inline T2 *GetMcTrack(T1 *reco)
   {
      if constexpr (std::is_same<T1, MpdTrack>::value) {
         if (reco->GetID() != 0) {
            return GetMCTrackDst(reco->GetID());
         } else
            return nullptr;
      } else {
         if (reco->mcTrackIndex() != 0) {
            return GetMCtrackMiniDst(reco->mcTrackIndex());
         }
         return nullptr;
      }
   }

   template <class T1, class T2>
   void MatchByDaughter(Int_t dau);
   template <class T1, class T2>
   void                   MatchByDaughters();
   inline MpdMCTrack *    GetMCTrackDst(Int_t index) const { return (MpdMCTrack *)fMcTracks->UncheckedAt(index); };
   inline MpdMiniMcTrack *GetMCtrackMiniDst(Int_t index) const
   {
      return (MpdMiniMcTrack *)fMiniMcTracks->UncheckedAt(index);
   }
   void               MatchV0ByMomentum();
   void               MatchByMc();
   void               MatchByMomentum();
   virtual InitStatus Init();

public:
   MpdV0Matcher(MpdV0::EParticleType pidHypo = MpdV0::EParticleType::kPdgHypo,
                EMatchType           type    = EMatchType::kMatchByBothDaughters);
   void SetPidHypo(MpdV0::EParticleType pidHypo) { fPidHipo = pidHypo; };
   void SetMatchingMethod(EMatchType type) { fMethod = type; };
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
