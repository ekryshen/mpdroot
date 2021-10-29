/*
 * MpdHbtTrack.h
 *
 *  Created on: 1 paź 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef NICAMPDHBTTRACK_H_
#define NICAMPDHBTTRACK_H_

#include "NicaExpTrack.h"
#include "NicaV0Track.h"

class NicaTrackTpcPads;

class NicaMpdHbtTrack : public NicaExpTrackHelix {
  NicaTrackTpcPads* fPads;

public:
  NicaMpdHbtTrack();
  NicaMpdHbtTrack(const NicaMpdHbtTrack& other);
  NicaTrackTpcPads* GetPadsInfo() const { return fPads; };
  NicaMpdHbtTrack& operator=(const NicaMpdHbtTrack& other);
  virtual void CopyData(NicaTrack* track);
  virtual ~NicaMpdHbtTrack();
  ClassDef(NicaMpdHbtTrack, 1)
};

#endif /* MPDROOT_NICA_MPD_FORMAT_HBTFORMAT_NICAMPDHBTTRACK_H_ */
