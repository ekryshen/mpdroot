//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFwdHit
/// \brief Reconstructed hit, reusing FairHit data members
/// \author Evgeny Kryshen (PNPI, Gatchina)
//------------------------------------------------------------------------------------------------------------------------

#ifndef __MPD_FWD_HIT_H
#define __MPD_FWD_HIT_H

#include "FairHit.h"

class MpdFwdHit : public FairHit {
public:
   MpdFwdHit();
   virtual ~MpdFwdHit();
   ClassDef(MpdFwdHit, 5)
};
#endif // #ifndef __MPD_FWD_HIT_H
