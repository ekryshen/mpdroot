/*
 * MpdV0CandicateCut.h
 *
 *  Created on: 16 gru 2021
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_PHYSICS_COMMON_V0_MPDV0CANDIDATECUT_H_
#define MPDROOT_PHYSICS_COMMON_V0_MPDV0CANDIDATECUT_H_

#include <Rtypes.h>
#include <RtypesCore.h>
#include <TObject.h>

class MpdV0Track;

class MpdV0CandidateCut : public TObject {
   Double_t fDca1to2;
   Double_t fDecayLenght;
   Double_t fDcaPrimVtx;

public:
   MpdV0CandidateCut() : fDca1to2(2.0), fDecayLenght(0), fDcaPrimVtx(1){};
   void SetMaxDcaDaughters(Double_t val) { fDca1to2 = val; }
   void SetMinDecayLenght(Double_t val) { fDecayLenght = val; }
   void SetMaxDcaPrimVtx(Double_t val) { fDcaPrimVtx = val; };
   virtual ~MpdV0CandidateCut(){};
   virtual Bool_t Pass(MpdV0Track &track);
   MpdV0CandidateCut(const MpdV0CandidateCut &other) = default;
   MpdV0CandidateCut &operator=(const MpdV0CandidateCut &other) = default;
   ClassDef(MpdV0CandidateCut, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_MPDV0CANDIDATECUT_H_ */
