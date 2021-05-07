/*
 * NicaMpdTrackTpcPads.h
 *
 *  Created on: 2 lip 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_FORMAT_MPDEVENT_NICAMPDTRACKTPCPADS_H_
#define MPDROOT_NICA_MPD_FORMAT_MPDEVENT_NICAMPDTRACKTPCPADS_H_

#include "NicaMpdTrack.h"
#include "NicaTpcSectorGeo.h"
#include "NicaTrackTpcPads.h"

class NicaMpdTrackTpcPads : public NicaMpdTrack {
 protected:
  NicaTrackTpcPads *fPads;
  /**
   *
   * @return true if pads were calculated
   */
  Bool_t PadsCalculated() const;

 public:
  NicaMpdTrackTpcPads();
  NicaMpdTrackTpcPads(const NicaMpdTrackTpcPads &other);
  NicaMpdTrackTpcPads &operator=(const NicaMpdTrackTpcPads &other);

  /**
   * calculate pads and pad id's this is usually used for pair-like cuts
   * @param shift - shift  by - vertex
   */
  void CalculatePads(Bool_t shift = kFALSE);
  NicaTrackTpcPads *GetTpcPadsInfo() const { return fPads; };
  virtual void Update(MpdTrack *track);
  virtual void CopyData(NicaTrack *other);
  virtual ~NicaMpdTrackTpcPads();
  ClassDef(NicaMpdTrackTpcPads, 1)
};

#endif /* MPDROOT_NICA_MPD_FORMAT_MPDEVENT_NICAMPDTRACKTPCPADS_H_ */
