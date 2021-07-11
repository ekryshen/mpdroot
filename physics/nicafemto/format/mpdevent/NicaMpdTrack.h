/*
 * NicaMpdTrack.h
 *
 *  Created on: 28 mar 2017
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef NICAMPDTRACK_H_
#define NICAMPDTRACK_H_

#include "MpdTrack.h"
#include "NicaExpTrack.h"
#include "NicaMpdTrack.h"
#include "NicaTrack.h"


class NicaMpdTrack : public NicaExpTrackHelix {
protected:
  TVector3 *fFirstPoint, *fLastPoint;
  NicaTpcTrack* fTpcTrack;
  NicaToFTrack* fToFTrack;
  ULong64_t fHitsMap;
  ULong64_t fSharedHitsMap;

public:
  NicaMpdTrack();
  NicaMpdTrack(const NicaMpdTrack& other);
  NicaMpdTrack& operator=(const NicaMpdTrack& other);
  inline TVector3* GetLastPoint() const { return fLastPoint; };
  inline TVector3* GetFirstPoint() const { return fFirstPoint; };
  inline NicaTpcTrack* GetTpcTrack() const { return fTpcTrack; };
  inline NicaToFTrack* GetToFTrack() const { return fToFTrack; };
  NicaDetectorTrack* GetDetTrack(const UInt_t detID) const;
  inline ULong64_t GetHitMap() const { return fHitsMap; };
  inline ULong64_t GetSharedHitMap() const { return fSharedHitsMap; };
  inline Bool_t HaveHit(Int_t layer) const { return TESTBIT(fHitsMap, layer); };
  inline Bool_t HaveSharedHit(Int_t layer) const { return TESTBIT(fSharedHitsMap, layer); };
  inline Int_t GetNHitsTpc() const { return fTpcTrack->GetNHits(); };
  virtual void Update(MpdTrack* track);
  void SetHitMap(ULong64_t hitMap) { fHitsMap = hitMap; };
  void PrintHitMaps() const;
  virtual void CopyData(NicaTrack* other);
  virtual ~NicaMpdTrack();
  ClassDef(NicaMpdTrack, 1)
};

#endif /* NICAMPDTRACK_H_ */
