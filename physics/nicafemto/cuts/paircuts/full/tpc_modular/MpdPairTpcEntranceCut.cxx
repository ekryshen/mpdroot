/*
 * MpdTpcEntranceCut.cxx
 *
 *  Created on: 30 paź 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdPairTpcEntranceCut.h"
namespace MpdPadsFormat {
MpdPairTpcEntranceCut::MpdPairTpcEntranceCut()
    : MpdNominalTpcPairPadsDistanceCut(1) {
  SetUnitName("TpcNominalEntrance [cm]");
  SetMinMax(0, 1E+5);
}

Bool_t MpdPairTpcEntranceCut::Pass(NicaTwoTrack* pair) {
  InitPass(pair);
  if (fFirstCommonPad == 0) {
    SetValue(GetDistance(0));
  } else {
    SetValue(TMath::TwoPi() * 40);
  }
  return Validate();
}

MpdPairTpcEntranceCut::~MpdPairTpcEntranceCut() {}

MpdPairTpcEntranceCut2D::MpdPairTpcEntranceCut2D()
    : MpdNominalTpcPairPadsDistanceCut(2) {
  SetUnitName("TpcNominalEntrance XY [cm]", XY());
  SetUnitName("TpcNominalEntrance Z [cm]", Z());
  SetMinMax(0, 1E+5, XY());
  SetMinMax(0, 1E+5, Z());
}

Bool_t MpdPairTpcEntranceCut2D::Pass(NicaTwoTrack* pair) {
  InitPass(pair);
  if (fFirstCommonPad == 0) {
    TVector3 dist = GetDistance3D(0);
    SetValue(dist.Pt(), XY());
    SetValue(TMath::Abs(dist.Z()), Z());
  } else {
    SetValue(TMath::TwoPi() * 40, 0);
    SetValue(TMath::TwoPi() * 40), 1;
  }
  return Validate();
}

MpdPairTpcEntranceCut2D::~MpdPairTpcEntranceCut2D() {}
}  // namespace MpdPadsFormat
