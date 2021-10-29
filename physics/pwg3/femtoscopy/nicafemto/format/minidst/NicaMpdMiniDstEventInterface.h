/*
 * NicaMpdMiniDstEventInterface.h
 *
 *  Created on: 17 kwi 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef NICAMPDMINIDSTEVENTINTERFACE_H_
#define NICAMPDMINIDSTEVENTINTERFACE_H_

#include <TVector3.h>

#include "MpdMiniEvent.h"
#include "NicaEventInterface.h"
#include "NicaTrackClones.h"

class NicaMpdMiniDstEvent;

class NicaMpdMiniDstEventInterface : public NicaEventInterface {
  friend class NicaMpdMiniDstEvent;
  friend class NicaMpdMiniDstEventV0;
  NicaTrackClones* fEvent;
  NicaTrackClones* fTracks;
  NicaTrackClones* fTofInfo;
  NicaTrackClones* fEmcInfo;

protected:
  virtual void ConnectToTree();
  void Register(Bool_t write);
  MpdMiniEvent* GetMiniEvent() { return (MpdMiniEvent*) fEvent->UncheckedAt(0); }

public:
  NicaMpdMiniDstEventInterface();
  virtual void CopyData(NicaEventInterface* s);
  virtual void Compress(Int_t* map, Int_t map_size);
  virtual void CopyAndCompress(NicaEventInterface* s, Int_t* map, Int_t map_size);
  virtual void FillTrackInterface(NicaTrackInterface* track, Int_t index);
  virtual Int_t GetTotalTrackNo() const;
  virtual TObject* GetRawEventPointer() const { return fEvent; };
  virtual NicaTrackInterface* GetTrackInterface() const;
  virtual TObject* GetRawTrackPointer(Int_t index) const;
  /** GETTERS AND SETTERS **/
  virtual void SetRunInfoId(Int_t i);
  virtual Int_t GetRunInfoId() const;
  virtual void SetMagneticField(TVector3 mag);
  virtual TVector3 GetMagneticField() const;
  virtual TLorentzVector GetVertexError() const;
  virtual TLorentzVector GetVertex() const;
  virtual ~NicaMpdMiniDstEventInterface();
  ClassDef(NicaMpdMiniDstEventInterface, 1)
};

#endif /* NICAMPDMINIDSTEVENTINTERFACE_H_ */
