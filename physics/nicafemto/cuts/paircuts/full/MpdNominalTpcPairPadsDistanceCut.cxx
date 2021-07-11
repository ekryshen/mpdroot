/*
 * MpdTpcDistanceCut.cxx
 *
 *  Created on: 30 paź 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */

#include "MpdNominalTpcPairPadsDistanceCut.h"

#include "NicaEvent.h"
#include "NicaMpdTrackTpcPads.h"
#include "NicaTwoTrack.h"

namespace MpdPadsFormat {
  MpdNominalTpcPairPadsDistanceCut::MpdNominalTpcPairPadsDistanceCut(Int_t size) :
    MpdTpcPadsPairCut(size),
    fCosDipAngle1(0),
    fCosDipAngle2(0),
    fSinDipAngle1(0),
    fSinDipAngle2(0),
    fCosPhase1(0),
    fCosPhase2(0),
    fSinPhase1(0),
    fSinPhase2(0),
    fTrack1(NULL),
    fTrack2(NULL),
    fdX(0),
    fdY(0),
    fdZ(0),
    fFirstCommonPad(0),
    fLastCommonPad(0) {}

  Double_t MpdNominalTpcPairPadsDistanceCut::GetDistance(Int_t lay) const {
    Double_t s1             = fTrack1->GetTpcPadsInfo()->GetPath(lay);
    const NicaHelix& helix1 = fTrack1->GetHelix();
    /* manual calcuations instead of helix, speed */
    Double_t x1 =
      helix1.GetStartX()
      + (TMath::Cos(helix1.GetPhi0() + s1 * helix1.GetH() * helix1.GetCurv() * fCosDipAngle1) - fCosPhase1) / helix1.GetCurv();
    Double_t y1 =
      helix1.GetStartY()
      + (TMath::Sin(helix1.GetPhi0() + s1 * helix1.GetH() * helix1.GetCurv() * fCosDipAngle1) - fSinPhase1) / helix1.GetCurv();
    Double_t z1             = helix1.GetStartZ() + s1 * fSinDipAngle1;
    Double_t s2             = fTrack2->GetTpcPadsInfo()->GetPath(lay);
    const NicaHelix& helix2 = fTrack2->GetHelix();
    Double_t x2 =
      helix2.GetStartX()
      + (TMath::Cos(helix2.GetPhi0() + s2 * helix2.GetH() * helix2.GetCurv() * fCosDipAngle2) - fCosPhase2) / helix2.GetCurv();
    Double_t y2 =
      helix2.GetStartY()
      + (TMath::Sin(helix2.GetPhi0() + s2 * helix2.GetH() * helix2.GetCurv() * fCosDipAngle2) - fSinPhase2) / helix2.GetCurv();
    Double_t z2 = helix2.GetStartZ() + s2 * fSinDipAngle2;
    x1          = x1 - fdX - x2;
    y1          = y1 - fdY - y2;
    z1          = z1 - fdZ - z2;
    return TMath::Sqrt(x1 * x1 + y1 * y1 + z1 * z1);
  }

  void MpdNominalTpcPairPadsDistanceCut::InitPass(NicaTwoTrack* pair) {
    fTrack1 = (NicaMpdTrackTpcPads*) pair->GetTrack1();
    fTrack2 = (NicaMpdTrackTpcPads*) pair->GetTrack2();
    fTrack1->CalculatePads();
    fTrack2->CalculatePads();
    fCosDipAngle1          = TMath::Cos(fTrack1->GetHelix().GetDipAngle());
    fCosDipAngle2          = TMath::Cos(fTrack2->GetHelix().GetDipAngle());
    fSinDipAngle1          = TMath::Sin(fTrack1->GetHelix().GetDipAngle());
    fSinDipAngle2          = TMath::Sin(fTrack2->GetHelix().GetDipAngle());
    fCosPhase1             = TMath::Cos(fTrack1->GetHelix().GetPhi0());
    fCosPhase2             = TMath::Cos(fTrack2->GetHelix().GetPhi0());
    fSinPhase1             = TMath::Sin(fTrack1->GetHelix().GetPhi0());
    fSinPhase2             = TMath::Sin(fTrack2->GetHelix().GetPhi0());
    fdX                    = fTrack1->GetEvent()->GetVertex()->X() - fTrack2->GetEvent()->GetVertex()->X();
    fdY                    = fTrack1->GetEvent()->GetVertex()->Y() - fTrack2->GetEvent()->GetVertex()->Y();
    fdZ                    = fTrack1->GetEvent()->GetVertex()->Z() - fTrack2->GetEvent()->GetVertex()->Z();
    const Int_t start1_pad = fTrack1->GetTpcPadsInfo()->GetFirstGoodPad();
    const Int_t start2_pad = fTrack2->GetTpcPadsInfo()->GetFirstGoodPad();
    const Int_t end1_pad   = fTrack1->GetTpcPadsInfo()->GetLastGoodPad() + 1;
    const Int_t end2_pad   = fTrack2->GetTpcPadsInfo()->GetLastGoodPad() + 1;
    fFirstCommonPad        = TMath::Max(start1_pad, start2_pad);
    fLastCommonPad         = TMath::Min(end1_pad, end2_pad);
  }

  MpdNominalTpcPairPadsDistanceCut::MpdNominalTpcPairPadsDistanceCut(const MpdNominalTpcPairPadsDistanceCut& other) :
    MpdTpcPadsPairCut(other),
    fCosDipAngle1(other.fCosDipAngle1),
    fCosDipAngle2(other.fCosDipAngle2),
    fSinDipAngle1(other.fSinDipAngle1),
    fSinDipAngle2(other.fSinDipAngle2),
    fCosPhase1(other.fCosPhase1),
    fCosPhase2(other.fCosPhase2),
    fSinPhase1(other.fSinPhase1),
    fSinPhase2(other.fSinPhase2),
    fTrack1(other.fTrack1),
    fTrack2(other.fTrack2),
    fdX(other.fdX),
    fdY(other.fdY),
    fdZ(other.fdZ),
    fFirstCommonPad(other.fFirstCommonPad),
    fLastCommonPad(other.fLastCommonPad) {}

  MpdNominalTpcPairPadsDistanceCut& MpdNominalTpcPairPadsDistanceCut::operator=(const MpdNominalTpcPairPadsDistanceCut& other) {
    if (this != &other) {
      MpdTpcPadsPairCut::operator=(other);
      fCosDipAngle1              = other.fCosDipAngle1;
      fCosDipAngle2              = other.fCosDipAngle2;
      fSinDipAngle1              = other.fSinDipAngle1;
      fSinDipAngle2              = other.fSinDipAngle2;
      fCosPhase1                 = other.fCosPhase1;
      fCosPhase2                 = other.fCosPhase2;
      fSinPhase1                 = other.fSinPhase1;
      fSinPhase2                 = other.fSinPhase2;
      fFirstCommonPad            = other.fFirstCommonPad;
      fLastCommonPad             = other.fLastCommonPad;
      fTrack1                    = other.fTrack1;
      fTrack2                    = other.fTrack2;
      fdX                        = other.fdX;
      fdY                        = other.fdY;
      fdZ                        = other.fdZ;
    }
    return *this;
  }

  MpdNominalTpcPairPadsDistanceCut::~MpdNominalTpcPairPadsDistanceCut() {}

  TVector3 MpdNominalTpcPairPadsDistanceCut::GetDistance3D(Int_t lay) const {
    Double_t s1             = fTrack1->GetTpcPadsInfo()->GetPath(lay);
    const NicaHelix& helix1 = fTrack1->GetHelix();
    /* manual calcuations instead of helix, speed */
    Double_t x1 =
      helix1.GetStartX()
      + (TMath::Cos(helix1.GetPhi0() + s1 * helix1.GetH() * helix1.GetCurv() * fCosDipAngle1) - fCosPhase1) / helix1.GetCurv();
    Double_t y1 =
      helix1.GetStartY()
      + (TMath::Sin(helix1.GetPhi0() + s1 * helix1.GetH() * helix1.GetCurv() * fCosDipAngle1) - fSinPhase1) / helix1.GetCurv();
    Double_t z1             = helix1.GetStartZ() + s1 * fSinDipAngle1;
    Double_t s2             = fTrack2->GetTpcPadsInfo()->GetPath(lay);
    const NicaHelix& helix2 = fTrack2->GetHelix();
    Double_t x2 =
      helix2.GetStartX()
      + (TMath::Cos(helix2.GetPhi0() + s2 * helix2.GetH() * helix2.GetCurv() * fCosDipAngle2) - fCosPhase2) / helix2.GetCurv();
    Double_t y2 =
      helix2.GetStartY()
      + (TMath::Sin(helix2.GetPhi0() + s2 * helix2.GetH() * helix2.GetCurv() * fCosDipAngle2) - fSinPhase2) / helix2.GetCurv();
    Double_t z2 = helix2.GetStartZ() + s2 * fSinDipAngle2;
    x1          = x1 - fdX - x2;
    y1          = y1 - fdY - y2;
    z1          = z1 - fdZ - z2;
    return TVector3(x1, y1, z1);
  }
}  // namespace MpdPadsFormat
