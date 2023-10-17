//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFwdPoint
///
/// \brief
/// \author Evgeny Kryshen (PNPI, Gatchina)
//------------------------------------------------------------------------------------------------------------------------

#ifndef __MPD_FWD_HIT_H
#define __MPD_FWD_HIT_H

#include "FairHit.h"

class MpdFwdHit : public FairHit {
private:
  float fTime;
  float fTrackLength;
public:
   MpdFwdHit();
   void SetTime(float time) { fTime = time; }
   float GetTime() { return fTime; }
   void SetLength(float length) { fTrackLength = length; }
   float GetLength() { return fTrackLength; }
   
   virtual ~MpdFwdHit();
   ClassDef(MpdFwdHit, 3)
};
#endif // #ifndef __MPD_FWD_HIT_H
