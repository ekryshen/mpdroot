/*
 * NicaTpcPads.h
 *
 *  Created on: 1 paź 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef NICATRACKTPCPADS_H_
#define NICATRACKTPCPADS_H_

#include "NicaHelix.h"
#include "NicaMpdConst.h"
#include <TVector.h>

/**
 * class that represent information about pads assigned to tracl
 */
class TLorentzVector;

class NicaTrackTpcPads : public TObject {
  Int_t fPadID[NicaMpdConst::TpcLayers];
  Float_t fPaths[NicaMpdConst::TpcLayers];
  Short_t fPadsNo[2];
  NicaHelix fNominalHelix;

public:
  NicaTrackTpcPads();
  /**
   *
   * @return number of pads
   */
  Int_t GetPadsNo() const { return fPadsNo[1] - fPadsNo[0]; };
  /**
   *
   * @return layer number with first pad
   */
  Int_t GetFirstGoodPad() const { return fPadsNo[0]; };
  /**
   *
   * @return layer number with last good pad
   */
  Int_t GetLastGoodPad() const { return fPadsNo[1] - 1; };
  /**
   *
   * @return maximal number of layers
   */
  Int_t GetMaxPadsNo() const { return NicaMpdConst::TpcLayers; };
  /**
   *
   * @param lay
   * @return ID of pad at given layer
   */
  Int_t GetPadID(Int_t lay) const { return fPadID[lay]; };
  /**
   *
   * @param lay
   * @return helix path at given layer
   */
  Double_t GetPath(Int_t lay) const { return fPaths[lay]; };
  /**
   *
   * @return true if properties are calculated
   */
  Bool_t AreCalculated() const;
  /**
   *
   * @return nominal helix (if Calculate() was called with event vertex)
   */
  const NicaHelix& GetNominalHelix() const { return fNominalHelix; };
  /**
   *
   * @param lay
   * @return momentum calculated at give layer by using helix
   */
  TVector3 GetMomentumAtLayer(Int_t lay) const { return fNominalHelix.EvalMom(fPaths[lay]); }
  /**
   *
   * @param lay
   * @return momentum calculated at givel layer by using helix
   */
  TVector3 GetPositionAtLayer(Int_t lay) const { return fNominalHelix.EvalPos(fPaths[lay]); };
  /**
   *
   * @param lay
   * @return radius that correspond to given layer
   */
  Double_t GetR(Int_t lay) const;
  /**
   * set status of this object as "not calculated"
   */
  void Reset() { fPadsNo[0] = -2; };
  /**
   * calculate the properties of object
   * @param helix helix that describes this track
   * @param vector vertex position, set vector==null to calculate real pads ID's
   * or event vertex to calculate nominal pads
   */
  void Calculate(const NicaHelix& helix, TLorentzVector* vector);
  virtual ~NicaTrackTpcPads();
  NicaTrackTpcPads(const NicaTrackTpcPads& other);
  NicaTrackTpcPads& operator=(const NicaTrackTpcPads& other);
  ClassDef(NicaTrackTpcPads, 1)
};

#endif /* MPDROOT_NICA_MPD_FORMAT_MPDEVENT_NICATRACKTPCPADS_H_ */
