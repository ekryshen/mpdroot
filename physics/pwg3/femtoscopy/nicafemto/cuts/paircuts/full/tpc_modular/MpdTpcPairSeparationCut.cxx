/*
 * MpdTpcSeparationCut.cxx
 *
 *  Created on: 30 paź 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdTpcPairSeparationCut.h"
#include "NicaTwoTrackCut.h"
namespace MpdPadsFormat {
MpdTpcPairSeparationCut::MpdTpcPairSeparationCut()
    : MpdNominalTpcPairPadsDistanceCut(5) {
  SetUnitName("TpcEntrance [cm]", TpcEntrance());
  SetUnitName("TpcExit [cm]", TpcExit());
  SetUnitName("TpcAverage [cm]", TpcAverage());
  SetUnitName("TpcMinimal[cm", TpcMinimal());
  SetUnitName("TpcMaximal[cm]", TpcMaximal());
  for (int i = 0; i < 5; i++) {
    SetMinMax(0, 1E+6, i);
  }
}

Bool_t MpdTpcPairSeparationCut::Pass(NicaTwoTrack* pair) {
  InitPass(pair);
  Double_t average = 0;
  Double_t minimum = 1E+6;
  Double_t maximum = -1;
  SetValue(TMath::TwoPi() * 40, TpcEntrance());
  SetValue(TMath::TwoPi() * 140, TpcExit());
  for (int iLay = fFirstCommonPad; iLay < fLastCommonPad; iLay++) {
    Double_t dst = GetDistance(iLay);
    if (iLay == 0) SetValue(dst, TpcEntrance());
    if (iLay == 52) SetValue(dst, TpcExit());
    average += dst;
    minimum = TMath::Min(minimum, dst);
    maximum = TMath::Max(maximum, dst);
  }
  SetValue(minimum, TpcMinimal());
  SetValue(maximum, TpcMaximal());
  Double_t pads = fLastCommonPad - fFirstCommonPad;
  if (pads != 0)
    SetValue(average / pads, TpcAverage());
  else
    SetValue(1000, TpcAverage());
  return Validate();
}

MpdTpcPairSeparationCut::~MpdTpcPairSeparationCut() {}
}  // namespace MpdPadsFormat
