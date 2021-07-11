/*
 * NicaMpdMiniDstEventV0Interface.h
 *
 *  Created on: 25 lut 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_FORMAT_MINIDSTV0_NICAMPDMINIDSTEVENTV0INTERFACE_H_
#define MPDROOT_NICA_MPD_FORMAT_MINIDSTV0_NICAMPDMINIDSTEVENTV0INTERFACE_H_

#include "NicaMpdMiniDstEventInterface.h"

class NicaMpdMiniDstEventV0Interface : public NicaMpdMiniDstEventInterface {
  friend class NicaMpdMiniDstEventV0;
  NicaTrackClones* fV0Tracks;
  NicaTrackClones* fV0Links;

protected:
  virtual void ConnectToTree();
  void Register(Bool_t write);

public:
  NicaMpdMiniDstEventV0Interface();
  virtual void CopyData(NicaEventInterface* s);
  virtual void Compress(Int_t* map, Int_t map_size);
  virtual void CopyAndCompress(NicaEventInterface* s, Int_t* map, Int_t map_size);
  virtual NicaTrackInterface* GetTrackInterface() const;
  virtual Int_t GetTotalTrackNo() const;
  virtual TClonesArray* GetV0Links() const { return fV0Links->GetArray(); };
  virtual ~NicaMpdMiniDstEventV0Interface();
  ClassDef(NicaMpdMiniDstEventV0Interface, 1)
};


#endif /* MPDROOT_NICA_MPD_FORMAT_MINIDSTV0_NICAMPDMINIDSTEVENTV0INTERFACE_H_ */
