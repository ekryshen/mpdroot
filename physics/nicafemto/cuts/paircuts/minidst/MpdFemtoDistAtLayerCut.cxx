/*
 * MpdFemtoDistAtLayerCut.cxx
 *
 *  Created on: 21 lip 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */

#include "MpdFemtoDistAtLayerCut.h"

#include "NicaMpdConst.h"
#include "NicaMpdHbtTrack.h"
#include "NicaPackage.h"
#include "NicaParameter.h"
#include "NicaTrackTpcPads.h"
#include "NicaTwoTrack.h"

namespace MpdHbtDst {

  MpdFemtoDistAtLayerCut::MpdFemtoDistAtLayerCut() : MpdFemtoPairCut(3), fLayerNo(0) {
    SetUnitName("#DeltaXY [cm[", 0);
    SetUnitName("#deltaZ [cm]", 1);
    SetUnitName("#deltaXYZ [cm]", 2);
    SetMinMax(1, -1, 0);
    SetMinMax(1, -1, 1);
    SetMinMax(1, -1, 2);
  }

  void MpdFemtoDistAtLayerCut::SetLayerNo(Int_t lay) {
    if (lay >= 0 && lay < NicaMpdConst::TpcLayers) { fLayerNo = lay; }
  }

  Bool_t MpdFemtoDistAtLayerCut::Pass(NicaTwoTrack* pair) {
    NicaMpdHbtTrack* tr1    = (NicaMpdHbtTrack*) pair->GetTrack1();
    NicaMpdHbtTrack* tr2    = (NicaMpdHbtTrack*) pair->GetTrack2();
    NicaTrackTpcPads* pads1 = tr1->GetPadsInfo();
    NicaTrackTpcPads* pads2 = tr2->GetPadsInfo();
    if (fLayerNo > pads1->GetLastGoodPad() || fLayerNo > pads2->GetLastGoodPad()) {
      SetValue(1E+3, 0);
      SetValue(1E+3, 1);
      SetValue(1E+3, 2);
    } else {
      TVector3 pos1 = pads1->GetPositionAtLayer(fLayerNo);
      TVector3 pos2 = pads2->GetPositionAtLayer(fLayerNo);
      pos1          = pos1 - pos2;
      SetValue(pos1.Pt(), 0);
      SetValue(pos1.Z(), 1);
      SetValue(pos1.Mag(), 2);
    }
    return AntiValidate();
  }

  NicaPackage* MpdFemtoDistAtLayerCut::Report() const {
    NicaPackage* p = MpdFemtoPairCut::Report();
    p->AddObject(new NicaParameterInt("Layer No", fLayerNo));
    return p;
  }

  MpdFemtoDistAtLayerCut::~MpdFemtoDistAtLayerCut() {}

} /* namespace MpdHbtDst */
