/// \class NicaTpcSectorGeo
///
/// MPD TPC sector geometry description
/// \author Alexander Zinchenko (LHEP, JINR, Dubna)

#include "NicaTpcSectorGeo.h"

#include <RtypesCore.h>
#include <TMath.h>
#include <TMathBase.h>
#include <TString.h>
#include <cstdlib>
#include <iostream>

#include "NicaHelixBase.h"

using std::cout;
using std::endl;

NicaTpcSectorGeo* NicaTpcSectorGeo::fgTpcSec = 0x0;

//__________________________________________________________________________
NicaTpcSectorGeo* NicaTpcSectorGeo::Instance(TString geo) {
  /// Get pointer to the TPC sector singleton object
  if (!fgTpcSec) {
    fgTpcSec = new NicaTpcSectorGeo(geo);

    // Automatic deletion
    std::atexit(DestroyInstance);
  }
  return fgTpcSec;
}

//__________________________________________________________________________
void NicaTpcSectorGeo::Init(TString geo) {
  /// Initialization of sector parameters
  fDphi = TMath::TwoPi() / fgkNsect;
  fPhi0 = -fDphi / 2.;

  cout << "\n !!!!! ***** ***** ***** !!!!!" << endl;
  // if (sensVol->getShape() == "TUBE") {
  if (geo == "simple") {
    // Old sensitive volume
    /*
fNrows[0] = fNrows[1] = 50;

Double_t rMin = ((TGeoTube*)shape)->GetRmin();
Double_t rMax = ((TGeoTube*)shape)->GetRmax();
cout << " ***** TPC sensitive volume: " << sv->GetName() << " " << shape->ClassName() << " " << rMin << " " << rMax << " " <<
((TGeoTube*)shape)->GetDZ() << endl;

fYsec[0] = 36.; // for consistency with earlier results
fYsec[1] = fYsec[2] = 96.; // for consistency with earlier results
fPadH[0] = fPadH[1] = (fYsec[1] - fYsec[0]) / fNrows[0];
*/
  } else if (geo == "full") {
    // New sensitive volume
    // fYsec[0] = params->At(5);
    // fYsec[1] = params->At(6);
    fYsec[0] = 40.3;
    // fYsec[2] = ((TGeoPgon*)shape)->Rmax(0);
    fPadH[0]       = 1.2;  // approx. pad height in inner ROC region
    fPadH[1]       = 1.8;  // 1.2; // approx. pad height in outer ROC region
    fNrows[0]      = 27;   // 33;
    fNrows[1]      = 26;   // 33;
    Double_t dy    = fPadH[0] * fNrows[0] + fPadH[1] * fNrows[1];
    fYsec[2]       = fYsec[0] + dy;
    Double_t scale = (fYsec[2] - fYsec[0]) / dy;
    fPadH[0] *= scale;
    fPadH[1] *= scale;
    fYsec[1] = fYsec[0] + fPadH[0] * fNrows[0];
    fPhi0    = 15.0 / TMath::RadToDeg();
    // Extract sensitive volume Z-coordinates
    fZmin = 0.4;
    // fZmax = fZmin + 2 * ((TGeoTube*)shape)->GetDZ();
    fZmax = 170.;
    fPhi0 = -0.261799;
  }

  // fPadH = (fYsec[1] - fYsec[0]) / fNrows;
  cout << " ***** TPC sector params - inner radius: " << fYsec[0] << ", outer radius: " << fYsec[2] << ", boundary: " << fYsec[1]
       << endl;
  cout << " Number of sectors: " << fgkNsect << ", phi0: " << fPhi0 * TMath::RadToDeg() << ", numbers of padrows: " << fNrows[0]
       << " " << fNrows[1] << ", pad heights: " << fPadH[0] << " " << fPadH[1] << endl;
  std::cout << "Raw phi" << fPhi0 << " DPhi " << fDphi << std::endl;
  cout << " !!!!! ***** ***** ***** !!!!!" << endl;

  fPadW[0] = fPadW[1] = 0.5;  // pad widths
  // Numbers of pads in rows
  // Double_t tan = TMath::Tan(fDphi/2), dead = 1.35; // dead area on one side
  Double_t tan = TMath::Tan(fDphi / 2), dead = 0.15;  // dead area on one side
  fNPadsInRows = new Int_t[NofRows()];
  for (Int_t j = 0; j < fNrows[0]; ++j)
    fNPadsInRows[j] = Int_t((tan * (fYsec[0] + j * fPadH[0]) - dead) / fPadW[0]);
  for (Int_t j = 0; j < fNrows[1]; ++j)
    fNPadsInRows[j + fNrows[0]] = Int_t((tan * (fYsec[1] + j * fPadH[1]) - dead) / fPadW[1]);
}

//__________________________________________________________________________
Int_t NicaTpcSectorGeo::Global2Local(const Double_t* xyzGlob, Double_t* xyzLoc, Int_t iSec0) {
  /// Transform global coordinates to local (sector)

  Double_t safety = 0.01;
  xyzLoc[2]       = xyzGlob[2];
  Int_t iSec      = iSec0;
  if (iSec0 < 0) {
    // Find sector No.
    Double_t phGlob = TMath::ATan2(xyzGlob[1], xyzGlob[0]);
    if (phGlob < 0) phGlob += TMath::TwoPi();
    iSec = Int_t((phGlob - fPhi0) / fDphi);
    if (iSec == fgkNsect) iSec = 0;
  } else
    iSec = (iSec0 >> kSectorS) & kSectorM;
  Double_t phSec = iSec * fDphi;
  Double_t cosPh = TMath::Cos(phSec);
  Double_t sinPh = TMath::Sin(phSec);
  Double_t x     = xyzGlob[0] * cosPh + xyzGlob[1] * sinPh;
  Double_t y     = -xyzGlob[0] * sinPh + xyzGlob[1] * cosPh;
  /*
if (x < fYsec[0] + safety || x > fYsec[2] - safety){
  std::cout<<"TRANSFORM"<<std::endl;
  std::cout<<TMath::Sqrt(xyzGlob[0]*xyzGlob[0]+xyzGlob[1]*xyzGlob[1])<<std::endl;
  std::cout<<"\t"<<xyzGlob[0]<<" "<<xyzGlob[1]<<std::endl;
  std::cout<<x-fYsec[0]<<" "<<y<<std::endl;
  //std::cout<<"\t"<<x<<" "<<fYsec[0]<<" "<<fYsec[2]<<std::endl;
}*/
  if (x < fYsec[0] + safety || x > fYsec[2] - safety) return -1;  // outside sector in Y
  xyzLoc[0] = y;
  xyzLoc[1] = x - fYsec[0];
  // if (x < fYsec[0] + safety || x > fYsec[1] - safety) return -1; // outside sector in Y
  Int_t row = fNrows[0];
  if (x <= fYsec[1])
    row = Int_t(xyzLoc[1] / fPadH[0]);
  else
    row += Int_t((x - fYsec[1]) / fPadH[1]);
  return (iSec |= (row << kPadrowS));
}

//__________________________________________________________________________
Int_t NicaTpcSectorGeo::Global2Local(const TVector3& xyzGlob, TVector3& xyzLoc, Int_t iSec0) {
  /// Transform global coordinates to local (sector)

  Double_t xyz1[3], xyz2[3];

  xyzGlob.GetXYZ(xyz1);
  Int_t iSec = Global2Local(xyz1, xyz2, iSec0);
  xyzLoc.SetXYZ(xyz2[0], xyz2[1], xyz2[2]);
  return iSec;
}

//__________________________________________________________________________
void NicaTpcSectorGeo::Local2Global(Int_t iSec, const Double_t* xyzLoc, Double_t* xyzGlob) {
  /// Transform local coordinates of sector iSec to global

  xyzGlob[2]     = xyzLoc[2];
  Double_t phSec = iSec * fDphi;
  Double_t cosPh = TMath::Cos(phSec);
  Double_t sinPh = TMath::Sin(phSec);
  Double_t x     = xyzLoc[1] + fYsec[0];
  Double_t y     = xyzLoc[0];
  xyzGlob[0]     = x * cosPh - y * sinPh;
  xyzGlob[1]     = x * sinPh + y * cosPh;
}

//__________________________________________________________________________
void NicaTpcSectorGeo::Local2Global(Int_t iSec, const TVector3& xyzLoc, TVector3& xyzGlob) {
  /// Transform local coordinates of sector iSec to global

  Double_t xyz1[3], xyz2[3];

  xyzLoc.GetXYZ(xyz1);
  Local2Global(iSec, xyz1, xyz2);
  xyzGlob.SetXYZ(xyz2[0], xyz2[1], xyz2[2]);
}

//__________________________________________________________________________
TVector2 NicaTpcSectorGeo::LocalPadPosition(Int_t padID) {
  /// Return local pad position for padID

  Int_t row  = PadRow(padID);
  Double_t x = 0.0, y = 0.0;
  if (row < fNrows[0])
    y = fPadH[0] * (row + 0.5);
  else
    y = fYsec[1] - fYsec[0] + fPadH[1] * (row - fNrows[0] + 0.5);
  return TVector2(x, y);
}

//__________________________________________________________________________

void NicaTpcSectorGeo::PadID(Float_t xloc, Float_t yloc, UInt_t& row, UInt_t& pad, Float_t& yNext) {
  /// Compute row, pad and distance to the nearest row (+ or -)

  row                = fNrows[0];
  Double_t lowHeight = fYsec[1] - fYsec[0];
  if (yloc <= lowHeight) {
    row   = Int_t(yloc / fPadH[0]);
    yNext = (row + 1) * fPadH[0] - yloc + 0.000001;
    if (yNext > fPadH[0] / 2) yNext -= (fPadH[0] + 0.000002);
    pad = Int_t(xloc / fPadW[0] + fNPadsInRows[row]);
  } else {
    Int_t dRow = Int_t((yloc - lowHeight) / fPadH[1]);
    row += dRow;
    yNext = (dRow + 1) * fPadH[1] + lowHeight - yloc + 0.000001;
    if (yNext > fPadH[1] / 2) yNext -= (fPadH[1] + 0.000002);
    pad = Int_t(xloc / fPadW[1] + fNPadsInRows[row]);
  }
}

void NicaTpcSectorGeo::CalculatePads(const NicaHelix& helix, Float_t* Paths, Short_t* Pads, Bool_t shift) {
  fCHCosDipAngle = TMath::Cos(helix.GetDipAngle()) * helix.GetCurv() * helix.GetH();
  fSecCH         = 1.0 / fCHCosDipAngle;
  fSinDipAngle   = TMath::Sin(helix.GetDipAngle());
  fSinPhase      = TMath::Sin(helix.GetPhi0());
  fCosPhase      = TMath::Cos(helix.GetPhi0());
  fHelixR        = 1.0 / helix.GetCurv();
  fPeriod        = helix.GetPeriod();
  if (fPeriod == 0 || TMath::IsNaN(fPeriod)) {
    Pads[0] = -1;
    Pads[1] = -1;
    return;
  }
  for (int i = 0; i < 53; i++) {
    Paths[i] = -1;
  }
  Double_t s;
  Int_t sect      = 0;
  Int_t start_pad = FindFirstPad(helix, Paths, sect);
  if (start_pad == -1) {  // no first pad
    Pads[0] = -1;
    Pads[1] = -1;
    return;
  }
  Int_t end_pad = FindLastPad(helix, Paths, start_pad + 1, sect);
  Pads[0]       = start_pad;
  Pads[1]       = end_pad;
}

Bool_t NicaTpcSectorGeo::GetPad(const NicaHelix& helix, Double_t R, Double_t& S, Int_t& Sect) {
  S = 0;
  if (CheckSector(helix, R, S, Sect)) {
    return kTRUE;  // found in this sector
  }
  /*no found in this sector, we have to look into another sector, but basically
   helix might pass through many rows in many sector, lets find pad with shortest
   helix lenght by staring from neighbors pads
  */

  for (int i = 0; i < 6; i++) {
    Int_t sect_plus  = (Sect + i) % 12;
    Int_t sect_minus = (Sect - i) % 12;
    Double_t s_plus  = 0;
    Double_t s_minus = 0;
    Int_t opt        = 0;
    if (CheckSector(helix, R, s_plus, sect_plus)) opt += 1;
    if (CheckSector(helix, R, s_minus, sect_minus)) opt += 2;
    switch (opt) {
      case 0:  // no overlap continue
        break;
      case 1: {  // plus okay, minus not
        Sect = sect_plus;
        S    = s_plus;
        return kTRUE;
      } break;
      case 2: {  // minus okay, plus not
        Sect = sect_minus;
        S    = s_minus;
        return kTRUE;
      } break;
      case 3: {  // both okay we need to find minimum
        if (s_plus < s_minus) {
          Sect = sect_plus;
          S    = s_plus;
        } else {
          Sect = sect_minus;
          S    = s_minus;
        }
        return kTRUE;
      } break;
    }
  }
  // check opposite sector
  Sect += 6;
  if (CheckSector(helix, R, S, Sect)) { return kTRUE; }
  return kFALSE;
}

Bool_t NicaTpcSectorGeo::CheckSector(const NicaHelix& helix, Double_t R, Double_t& S, Int_t Sect) {
  Double_t mz = Sect;
  if (Sect >= 12) mz = Sect - 12;
  const Double_t rot = TVector2::Phi_mpi_pi(TMath::TwoPi() * mz / 12.0);
  const Double_t phi = TVector2::Phi_mpi_pi(helix.GetPhi0() + rot);
  const Double_t x   = helix.GetStartX() * TMath::Cos(-rot) + helix.GetStartY() * TMath::Sin(-rot);
  Double_t magfunc   = TMath::ACos(-helix.GetCurv() * x + helix.GetCurv() * R + TMath::Cos(phi));
  if (magfunc == 0) return kFALSE;  // cannot find overlaped lines
  Double_t s1 = fSecCH * (-magfunc - phi);
  Double_t s2 = fSecCH * (magfunc - phi);
  while (s1 < 0)
    s1 += fPeriod;
  while (s2 < 0)
    s2 += fPeriod;
  S = TMath::Min(s1, s2);
  if (S == 0) return kFALSE;
  //	if(S>fPeriod*0.5)
  //	return kFALSE;
  Double_t X      = helix.GetStartX() + fHelixR * (TMath::Cos(helix.GetPhi0() + S * fCHCosDipAngle) - fCosPhase);
  Double_t Y      = helix.GetStartY() + fHelixR * (TMath::Sin(helix.GetPhi0() + S * fCHCosDipAngle) - fSinPhase);
  Double_t newPhi = TVector2::Phi_mpi_pi(TMath::ATan2(Y, X) + rot);

  if (newPhi < -fSectAngle) return kFALSE;
  if (newPhi > fSectAngle) return kFALSE;

  /*TVector3 pos = helix->Evaluate(S);
  pos.RotateZ(rot);
  if(TMath::Abs(pos.X()-R)>1E-13){
          std::cout<<" BAD CALC "<<R<<" "<<pos.X()<<" \t"<<(pos.X()-R)<<std::endl;
          std::cout<<"S="<<S<<"\t sector: "<<Sect<<" p = "<<X<<" "<<Y<<std::endl;
          std::cout<<"\tstart="<<helix->GetStartX()<<" "<<helix->GetStartY()<<std::endl;
          //helix->Print();
  }*/

  Double_t Z = helix.GetStartZ() + S * fSinDipAngle;
  if (TMath::Abs(Z) > GetZmax()) return kFALSE;
  return kTRUE;
}


Int_t NicaTpcSectorGeo::FindFirstPad(const NicaHelix& helix, Float_t* Paths, Int_t& sect) {
  Double_t Rmin = GetRocY(0);
  Double_t s    = 0;
  Double_t S[fgkNsect];
  for (int i = 0; i < NofRowsReg(0); i++) {  // loop over inner layers
    Double_t R               = Rmin + (0.5 + (Double_t) i) * PadHeight(0);
    Bool_t found_good_sector = kFALSE;
    for (int j = 0; j < fgkNsect; j++) {  // check all sectors, find all overlaps
      if (!CheckSector(helix, R, S[j], j)) {
        S[j] = 1E+39;
      } else {
        if (AdditionalCheck(helix, S[j], R, j)) {
          found_good_sector = kTRUE;
        } else {
          S[j] = 1E+39;
        }
      }
    }
    if (found_good_sector) {  // at least one good secor found
      Int_t min_sec  = fgkNsect;
      Double_t s_min = 1E+39;
      for (int j = 0; j < fgkNsect; j++) {
        if (S[j] < s_min) {
          s_min   = S[j];
          min_sec = j;
        }
      }
      Paths[i] = s_min;
      sect     = min_sec;
      return i;
    }
  }
  Rmin = GetRocY(1);
  for (int i = NofRowsReg(0); i < NofRows(); i++) {  // loop over outer layers
    Double_t R               = (0.5 + (Double_t)(i - NofRowsReg(0))) * PadHeight(1) + Rmin;
    Bool_t found_good_sector = kFALSE;
    for (int j = 0; j < fgkNsect; j++) {  // check all sectors, find all overlaps
      if (!CheckSector(helix, R, S[j], j)) {
        S[j] = 1E+39;
      } else {
        if (AdditionalCheck(helix, S[j], R, j)) {
          found_good_sector = kTRUE;
        } else {
          S[j] = 1E+39;
        }
      }
    }
    if (found_good_sector) {  // at least one good secor found
      Int_t min_sec  = fgkNsect;
      Double_t s_min = 1E+39;
      for (int j = 0; j < fgkNsect; j++) {
        if (S[j] < s_min) {
          s_min   = S[j];
          min_sec = j;
        }
      }
      Paths[i] = s_min;
      sect     = min_sec;
      return i;
    }
  }
  return -1;
}

Bool_t NicaTpcSectorGeo::AdditionalCheck(const NicaHelix& helix, Double_t S, Double_t R, Double_t sect) const {
  Double_t angle = TMath::TwoPi() * sect / ((Double_t) fgkNsect);
  TVector3 pos   = helix.EvalPos(S);
  pos.RotateZ(angle);
  TVector3 start = helix.GetStartPoint();
  if (TMath::Abs(pos.X() - R) > 1E-6) { return kFALSE; }
  // reject point that pass pad but is closer than start vector
  if (start.Pt() > R) return kFALSE;
  return kTRUE;
}

Int_t NicaTpcSectorGeo::FindLastPad(const NicaHelix& helix, Float_t* Paths, Int_t first_pad, Int_t& sect) {
  Double_t Rmin = GetRocY(0);
  Double_t s    = 0;
  for (int i = first_pad; i < NofRowsReg(0); i++) {  // loop over inner layers
    Double_t R = Rmin + (0.5 + (Double_t) i) * PadHeight(0);
    if (GetPad(helix, R, s, sect)) {
      Paths[i] = s;
    } else {
      return i;
    }
  }
  Rmin = GetRocY(1);
  if (first_pad < NofRowsReg(0)) first_pad = NofRowsReg(0);
  for (int i = first_pad; i < NofRows(); i++) {  // loop over outer layers
    Double_t R = (0.5 + (Double_t)(i - NofRowsReg(0))) * PadHeight(1) + Rmin;
    if (GetPad(helix, R, s, sect)) {
      Paths[i] = s;
    } else {
      return i;
    }
  }
  return 53;
}

//__________________________________________________________________________
