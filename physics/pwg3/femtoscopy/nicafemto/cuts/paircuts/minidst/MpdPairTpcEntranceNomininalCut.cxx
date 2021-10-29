/*
 * MpdPairTpcEntranceNomininalCut.cxx
 *
 *  Created on: 1 paź 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdPairTpcEntranceNomininalCut.h"
#include "NicaMpdHbtTrack.h"
#include "NicaTrackTpcPads.h"
#include "NicaTwoTrack.h"

namespace MpdHbtDst {
  MpdPairTpcEntranceNomininalCut::MpdPairTpcEntranceNomininalCut() : MpdFemtoPairCut(3) {
    SetUnitName("Nominal TPC entrance sep [cm]", XYZ());
    SetUnitName("Nominal TPC entrance sep XY[cm]", XY());
    SetUnitName("Nominal TPC entrance sep Z [cm]", Z());
  }

  Bool_t MpdPairTpcEntranceNomininalCut::Pass(NicaTwoTrack* pair) {
    NicaMpdHbtTrack* track1 = (NicaMpdHbtTrack*) pair->GetTrack1();
    NicaMpdHbtTrack* track2 = (NicaMpdHbtTrack*) pair->GetTrack2();
    Double_t s1, s2;
    track1->GetPadsInfo()->GetNominalHelix().PathLength(NicaMpdConst::TpcInnerDriftRadius, s1, s2);
    Double_t S1 = TMath::Min(s1, s2);
    if (s1 < 0) S1 = s2;
    if (s2 < 0) S1 = s1;
    track2->GetPadsInfo()->GetNominalHelix().PathLength(NicaMpdConst::TpcInnerDriftRadius, s1, s2);
    Double_t S2 = TMath::Min(s1, s2);
    if (s1 < 0) S2 = s2;
    if (s2 < 0) S2 = s1;
    TVector3 pos1 = track1->GetPadsInfo()->GetNominalHelix().EvalPos(S1);
    TVector3 pos2 = track2->GetPadsInfo()->GetNominalHelix().EvalPos(S2);
    pos2          = pos2 - pos1;
    SetValue(pos2.Mag(), XYZ());
    SetValue(pos2.Pt(), XY());
    SetValue(TMath::Abs(pos2.Z()), Z());
    return Validate();
  }

  MpdPairTpcEntranceNomininalCut::~MpdPairTpcEntranceNomininalCut() {}

}  // namespace MpdHbtDst
