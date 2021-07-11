/*
 * NicaMpdV0Matcher.h
 *
 *  Created on: 27 lut 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_TASKS_NICAMPDV0MATCHER_H_
#define MPDROOT_NICA_MPD_TASKS_NICAMPDV0MATCHER_H_

#include "FairTask.h"

class NicaTrackClones;
class MpdMiniMcTrack;
class NicaV0Track;
/**
 * class for match v0 with mc tracks in minidst
 */

class NicaMpdV0Matcher : public FairTask {
protected:
  virtual InitStatus Init();
  NicaTrackClones* fV0Tracks;
  NicaTrackClones* fV0Links;
  NicaTrackClones* fMcTracks;
  NicaTrackClones* fReTracks;
  Int_t fLastPrimary;
  Int_t FindMotherIndex(MpdMiniMcTrack* mc, NicaV0Track* v0);
  Bool_t fWrite;

public:
  NicaMpdV0Matcher();
  void Write(Bool_t wtr = kTRUE) { fWrite = wtr; };
  virtual void Exec(Option_t* opt);
  virtual ~NicaMpdV0Matcher();
  ClassDef(NicaMpdV0Matcher, 1)
};


#endif /* MPDROOT_NICA_MPD_TASKS_NICAMPDV0MATCHER_H_ */
