/*
 * MpdV0FinderBasic.h
 *
 *  Created on: 16 gru 2021
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_PHYSICS_COMMON_V0_MPDV0FINDERBASIC_H_
#define MPDROOT_PHYSICS_COMMON_V0_MPDV0FINDERBASIC_H_

#include <FairTask.h>
#include <Rtypes.h>
#include <RtypesCore.h>
#include <memory>
#include <MpdHelix.h>

class MpdV0CandidateCut;
class MpdV0DaughterCut;

/**
 * abstract class for finding vo particles
 */

class MpdEvent;

class MpdV0FinderBasic : public FairTask {
protected:
   enum class EFormat { kDst, kMiniDst };
   Bool_t                                   fInit;
   Bool_t                                   fWrite;
   Bool_t                                   fFirstV0;
   EFormat                                  fFormat;
   MpdEvent *                               fMpdEvent;
   TClonesArray *                           fMiniEvents;
   TClonesArray *                           fMiniTracks;
   TClonesArray *                           fV0s;
   std::vector<std::pair<TObject *, Int_t>> fFirstDaughters;
   std::vector<std::pair<TObject *, Int_t>> fSecondDaughters;
   virtual void                             ExecDst(Option_t *option)     = 0;
   virtual void                             ExecMiniDst(Option_t *option) = 0;
   virtual InitStatus                       Init();
   MpdV0DaughterCut *                       fFirstDaughterCut;
   MpdV0DaughterCut *                       fSecondDaughterCut;
   MpdV0CandidateCut *                      fCandicateCut;

public:
   MpdV0FinderBasic();
   void         SaveV0s(Bool_t write) { fWrite = write; };
   virtual void Exec(Option_t *option);
   virtual ~MpdV0FinderBasic();
   void SetFirstDaughterCut(const MpdV0DaughterCut &cut);
   void SetSecondDaughterCut(const MpdV0DaughterCut &cut);
   void SetCandicateCut(const MpdV0CandidateCut &cut);
   MpdV0FinderBasic(const MpdV0FinderBasic &other);
   MpdV0FinderBasic &operator=(const MpdV0FinderBasic &other);
   ClassDef(MpdV0FinderBasic, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_MPDV0FINDERBASIC_H_ */
