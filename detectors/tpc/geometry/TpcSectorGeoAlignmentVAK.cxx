/// TpcSectorGeoAlignmentVAK source file
///
/// \author Alexander Zinchenko (LHEP, JINR, Dubna)
/// \corrected by Valentin Kuzmin, SINP of Moscow State University

#include "TpcSectorGeoAlignmentVAK.h"
#include <iostream>
#include "TpcGas.h"

#include <TGeoManager.h>
#include <TGeoPgon.h>
#include <TGeoTube.h>
#include <TMath.h>
#include <TSystem.h>
#include <Riostream.h>

using std::cout;
using std::endl;

//______________________________________________________________________

void TpcSectorGeoAlignmentVAK::Init()
{
   /// Initialization of sectors parameters

   fZmin     = 1;
   fZmax     = 168.;
   fPadH[0]  = 1.2; // pad height in inner ROC region
   fPadH[1]  = 1.8; // pad height in outer ROC region
   fNrows[0] = 27;
   fNrows[1] = 26;
   fFNRows   = fNrows[0] + fNrows[1];
   flYsec[0] = 0.;
   fYsec[0]  = 40.3;
   flYsec[1] = flYsec[0] + fPadH[0] * fNrows[0];
   fYsec[1]  = fYsec[0] + flYsec[1];
   flYsec[2] = flYsec[1] + fPadH[1] * fNrows[1];
   fYsec[2]  = fYsec[0] + flYsec[2];
   fPadW[0] = fPadW[1] = 0.5; // pad widths
   fPadZ[0]            = 170.;
   fPadZ[1]            = -170.;
   fPadS[0]            = fPadW[0] * fPadH[0];
   fPadS[1]            = fPadW[1] * fPadH[1];
   N_InnerPads_min     = 40; // Pads number in the lowest row
   Bool_t *addInnerPadsRow =
      new Bool_t[fNrows[0]]{0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1};

   fDphi = TMath::TwoPi() / fgkNsect;
   fPhi0 = -fDphi / 2.;

   fNPadsInRow  = new Int_t[fFNRows]{0}; // Number of pads in the row
   fNPadsInRows = new Int_t[fFNRows]{0}; // Number of pads in all rows below
                                         // this one
   Int_t N_PadsinRow = N_InnerPads_min;
   fNPadsInRows[0]   = 0;
   fNPadsInRow[0]    = N_InnerPads_min;
   for (UInt_t k = 1; k < fNrows[0]; ++k) {
      fNPadsInRows[k] = fNPadsInRows[k - 1] + fNPadsInRow[k - 1];
      fNPadsInRow[k] += (addInnerPadsRow[k]) ? fNPadsInRow[k - 1] + 2 : fNPadsInRow[k - 1];
   };
   for (UInt_t k = 0; k < fNrows[1]; ++k) {
      fNPadsInRows[fNrows[0] + k] = fNPadsInRows[fNrows[0] + k - 1] + fNPadsInRow[fNrows[0] + k - 1];
      fNPadsInRow[fNrows[0] + k]  = fNPadsInRow[fNrows[0] + k - 1] + 2;
   };
   fNPads   = fNPadsInRows[fFNRows - 1] + fNPadsInRow[fFNRows - 1];
   fXsec[0] = 0.5 * fPadW[0] * fNPadsInRow[0];
   fXsec[1] = 0.5 * fPadW[0] * fNPadsInRow[fFNRows - 1];

   // Gas parameters
   // std::string tpcGasFile = gSystem->Getenv("VMCWORKDIR");
   // tpcGasFile += "/geometry/Ar-90_CH4-10.asc";
   // fGas = new TpcGas(tpcGasFile, 130);
   // fTimeMax = fZmax / fGas->VDrift();
   fNTimeBins = 512; // max number of time bins
   // fTimeBin = fTimeMax / fNTimeBins;
   fTimeBin        = 100; // 100 ns
   Double_t VDrift = fZmax / (fNTimeBins * fTimeBin * 1.E-9);
   fTimeMax        = (fZmax - fZmin) / VDrift;
   fTimeBinMax     = fTimeMax / fTimeBin;
   fZ2TimeBin      = fTimeBinMax / fZmax;
}

//__________________________________________________________________

void TpcSectorGeoAlignmentVAK::PrintSectorParam()
{

   cout << "************** TPC sector parameters *************" << endl;
   cout << " Number of sectors: " << 2 * fgkNsect << "\n"
        << "phi0: " << fPhi0 * TMath::RadToDeg() << "\n"
        << "numbers of padrows: " << fNrows[0] << " " << fNrows[1] << "\n"
        << "pad heights: " << fPadH[0] << " & " << fPadH[1] << "\n"
        << "Number of inner pads: " << fNPadsInRows[fNrows[0]] << "\n"
        << "Number of outer pads: " << fNPads - fNPadsInRows[fNrows[0]] << "\n"
        << "Full number of pads in the sector: " << fNPads << "\n"
        << "Full number of pads: " << 2 * fgkNsect * fNPads << endl;

   cout << "----------  Inner Pads ----------"
        << "\n"
        << "row     number of pads " << endl;
   for (UInt_t k = 0; k < fFNRows; ++k) {
      if (k < 10) {
         cout << " " << k << "         " << fNPadsInRow[k] << "\n";
         continue;
      }
      if (k == 27) {
         cout << "----------  Outer Pads ----------" << endl;
      }
      if (k < 40) {
         cout << k << "         " << fNPadsInRow[k] << "\n";
      } else {
         cout << k << "        " << fNPadsInRow[k] << "\n";
      }
   };
   cout << "************ end TPC sector params ***********\n" << endl;
}

//__________________________________________________________________

void TpcSectorGeoAlignmentVAK::GetSectorAlignment(TString *inFile)
{
   // the current sectors allignment from MPDS DB
   // !!!<= lowHeight
   // it should be read but at the moment if true
   // we put all null values (the ideal case of the allignment)
   for (Int_t i = 0; i < 24; i++) {
      aR0_lsc[i].SetXYZ(0, 0, 0);
      alpha_E[i] = 0;
      beta_E[i]  = 0;
      gamma_E[i] = 0;
      A_l2tl[i].SetToIdentity();
      A_tl2l[i].SetToIdentity();
   }

   // !!!>
}

//______________________________________________________________________

/// Get the sector allignment as a shift and Euler's angles
void TpcSectorGeoAlignmentVAK::GetAlignment(Int_t iSec, TVector3 &r0_loc, Double_t &alpha, Double_t &beta,
                                            Double_t &gamma)
{
   r0_loc = aR0_lsc[iSec];
   alpha  = alpha_E[iSec];
   beta   = beta_E[iSec];
   gamma  = gamma_E[iSec];
}

//______________________________________________________________________

// void TpcSectorGeoAlignmentVAK::TpcSectorGeoA(Double_t r0_lsc[3][24], Double_t* alpha,
//                                        Double_t* beta,Double_t* gamma) {
void TpcSectorGeoAlignmentVAK::TpcSectorGeoA(Double_t r0_lsc[3][24], Double_t alpha[24], Double_t beta[24],
                                             Double_t gamma[24])
{
   // use a special allignment
   for (Int_t i = 0; i < 24; i++) {
      aR0_lsc[i] = TVector3(r0_lsc[0][i], r0_lsc[1][i], r0_lsc[2][i]);
      alpha_E[i] = alpha[i];
      beta_E[i]  = beta[i];
      gamma_E[i] = gamma[i];
   };
   PutSectorAlignment();
};

//______________________________________________________________________

void TpcSectorGeoAlignmentVAK::PutSectorAlignment()
{ // Calculation of the final transformation with the allignment
   for (Int_t i = 0; i < 24; i++) {
      A_l2tl[i].SetToIdentity(); // Identity transformation
      // abg
      A_l2tl[i].RotateX(gamma_E[i]);    // Turn around X
      A_l2tl[i].RotateY(beta_E[i]);     // Turn around Y
      A_l2tl[i].RotateZ(alpha_E[i]);    // Turn around Z final: LCS->A_l2tl[
      A_tl2l[i]  = A_l2tl[i].Inverse(); //   A_l2tl[->LCS
      aG0_lsc[i] = -(A_tl2l[i] * aR0_lsc[i]);
   }
   for (Int_t i = 0; i < 24; i++) {
      Int_t    ir;
      Double_t z0loc, phi, tangle;
      T_tl2glo[i].SetToIdentity(); // Identity transformation
      if (i < 12) {
         ir     = i;
         phi    = ir * fDphi;
         tangle = -TMath::Pi() / 2 + phi; // rotation angle from GCS->A_l2tl[ (z>0)
         z0loc  = fPadZ[0];
      } else {
         T_tl2glo[i].RotateY(TMath::Pi()); // Turn Z in opposite direction
         ir     = i - 12;
         phi    = TMath::Pi() - ir * fDphi; // rotation angle from GCS->A_l2tl[ (z<0)
         tangle = phi - TMath::Pi() / 2;
         z0loc  = fPadZ[1];
      };

      T_tl2glo[i].RotateZ(tangle); // Turn around Z to the sector i
      // coordinates of the center of LCS of the sector i
      R0_tlsc[i].SetXYZ(fYsec[0] * TMath::Cos(phi), fYsec[0] * TMath::Sin(phi), z0loc);
      T_glo2tl[i] = T_tl2glo[i].Inverse();
      G0_tlsc[i]  = -(T_glo2tl[i] * R0_tlsc[i]);
      // coordinates of the center of LCS of the sector i with alignment
      R0_lsc[i] = R0_tlsc[i] + A_l2tl[i] * aR0_lsc[i];
      // Rotation matrix GSC->IdealLCS->LCS
      R_loc2glo[i] = T_tl2glo[i] * A_l2tl[i];
      // Rotation matrix (GSC->LSC) with alignment
      R_glo2loc[i] = R_loc2glo[i].Inverse();
      // coordinates of the center of GCS in LCS of the sector i with alignment
      G0_gsc[i] = -(R_glo2loc[i] * R0_lsc[i]);
   };
}

//______________________________________________________________________

void TpcSectorGeoAlignmentVAK::PutSectorAlignment(Int_t iSec, TVector3 R0, Double_t alpha, Double_t beta,
                                                  Double_t gamma)
{ // Calculation of the final transformation with the non zero allignment
   // only for one sector
   aR0_lsc[iSec] = R0;
   alpha_E[iSec] = alpha;
   beta_E[iSec]  = beta;
   gamma_E[iSec] = gamma;
   A_l2tl[iSec].SetToIdentity(); // Identity transformation
   // (-)gamma->(-)beta->(-)alpha transformation (back 2
   A_l2tl[iSec].RotateZ(alpha); // Turn around Z
   A_l2tl[iSec].RotateY(beta);  // Turn around Y
   A_l2tl[iSec].RotateX(gamma); // Turn around X
   A_tl2l[iSec] = A_l2tl[iSec].Inverse();
   Int_t    ir;
   Double_t z0loc, phi, tangle;
   T_tl2glo[iSec].SetToIdentity(); // Identity transformation
   if (iSec < 12) {
      ir     = iSec;
      tangle = -TMath::Pi() / 2 + ir * fDphi;
      z0loc  = fPadZ[0];
   } else {
      T_tl2glo[iSec].RotateY(TMath::Pi()); // Turn Z in opposite direction
      ir     = iSec - 12;
      tangle = TMath::Pi() / 2 - ir * fDphi;
      z0loc  = fPadZ[1];
   }
   T_tl2glo[iSec].RotateZ(tangle); // Turn aroun Z to the sector i
   // coordinates of the center of LCS of the sector iSec
   R0_lsc[iSec].SetXYZ(fYsec[0] * TMath::Cos(phi), fYsec[0] * TMath::Sin(phi), z0loc);
   // coordinates of the center of LCS of the sector iSec with alignment
   R0_lsc[iSec] += T_tl2glo[iSec] * R0;
   // Rotation matrix (LSC->GSC) with alignment
   R_loc2glo[iSec] = T_tl2glo[iSec] * A_l2tl[iSec];
   // Rotation matrix (GSC->LSC) with alignment
   R_glo2loc[iSec] = R_loc2glo[iSec].Inverse();
   //  coordinates of the center of GCS in LCS of the sector iSec with alignment
   G0_gsc[iSec] = -(R_glo2loc[iSec] * R0_lsc[iSec]);
}

//______________________________________________________________________

void TpcSectorGeoAlignmentVAK::PrintSectorAlignment(Int_t isec)
{
   cout << "TpcSectorGeoAlignmentVAK: ========= Alignment ================" << endl;
   Int_t s1, s2;
   if (isec < 0) {
      s1 = 0;
      s2 = 24;
   } else if (isec < 24) {
      s1 = isec;
      s2 = isec + 1;
   } else {
      cout << "TpcSectorGeoAlignmentVAK: No " << isec << "-th sector" << endl;
      return;
   }
   for (Int_t i = s1; i < s2; i++) {
      cout << "============= sector=" << i << " =============" << endl;
      printf("Euler's angles (%14e%14e%14e)\n", alpha_E[i], beta_E[i], gamma_E[i]);
      //   LCS->GCS
      printf("- - LCS->GCS:  s0(%14e%14e%14e)\n", R0_lsc[i](0), R0_lsc[i](1), R0_lsc[i](2));
      for (Int_t k = 0; k < 3; k++)
         printf("%14e%14e%14e\n", R_loc2glo[i](k, 0), R_loc2glo[i](k, 1), R_loc2glo[i](k, 2));
      //   TLCS->GCS
      printf("- - TLCS->GCS: s0(%14e%14e%14e)\n", R0_lsc[i](0), R0_lsc[i](1), R0_lsc[i](2));
      for (Int_t k = 0; k < 3; k++) printf("%14e%14e%14e\n", T_tl2glo[i](k, 0), T_tl2glo[i](k, 1), T_tl2glo[i](k, 2));
      //   LCS->TLCS
      printf("- - LCS->TLCS: s0(%14e%14e%14e)\n", aR0_lsc[i](0), aR0_lsc[i](1), aR0_lsc[i](2));
      for (Int_t k = 0; k < 3; k++) printf("%14e%14e%14e\n", A_l2tl[i](k, 0), A_l2tl[i](k, 1), A_l2tl[i](k, 2));
   }
   cout << "  " << endl;
}

//______________________________________________________________________

void TpcSectorGeoAlignmentVAK::PrintAlignmentTransformation()
{ // Print transformation matrixes due to allignment
   cout << "************** TPC sectors allignment *************" << endl;
   for (Int_t i = 0; i < 24; i++) {
      cout << " -----  sector: " << i << " -----" << endl;
      cout << "shift: (" << aR0_lsc[i].X() << "," << aR0_lsc[i].Y() << "," << aR0_lsc[i].z() << ")" << endl;
      cout << "Euler angles: alpha, beta, gamma (deg)"
           << "       " << TMath::RadToDeg() * alpha_E[i] << endl
           << "       " << TMath::RadToDeg() * beta_E[i] << endl
           << "       " << TMath::RadToDeg() * gamma_E[i] << endl
           << endl;
      cout << "transformation : (" << A_l2tl[i].XX() << "," << A_l2tl[i].XY() << "," << A_l2tl[i].XZ() << ")\n";
      cout << "                 (" << A_l2tl[i].YX() << "," << A_l2tl[i].YY() << "," << A_l2tl[i].YZ() << ")\n";
      cout << "                 (" << A_l2tl[i].ZX() << "," << A_l2tl[i].ZY() << "," << A_l2tl[i].ZZ() << ")\n\n";
   }
}
//______________________________________________________________________

void TpcSectorGeoAlignmentVAK::SectorTransformation(const Int_t iSec, TVector3 &SecCoor_shift, TRotation &toGloSecCoor)
{
   SecCoor_shift = R0_lsc[iSec];
   toGloSecCoor  = R_loc2glo[iSec];
}

//______________________________________________________________________

void TpcSectorGeoAlignmentVAK::SectorBackTransformation(const Int_t iSec, TVector3 &SecCoor_shift_toloc,
                                                        TRotation &toLocSecCoor)
{
   SecCoor_shift_toloc = G0_gsc[iSec];
   toLocSecCoor        = R_glo2loc[iSec];
}

//______________________________________________________________________

void TpcSectorGeoAlignmentVAK::LCS2TLCStransformation(const Int_t iSec, TVector3 &LCS2TLCSshift, TRotation &LCS2TLCS)
{
   LCS2TLCSshift = aR0_lsc[iSec];
   LCS2TLCS      = A_l2tl[iSec];
}

//______________________________________________________________________

void TpcSectorGeoAlignmentVAK::TLCS2LCStransformation(const Int_t iSec, TVector3 &TLCS2LCSshift, TRotation &TLCS2LCS)
{
   TLCS2LCSshift = aG0_lsc[iSec];
   TLCS2LCS      = A_tl2l[iSec];
}

//______________________________________________________________________

void TpcSectorGeoAlignmentVAK::TLCS2GCStransformation(const Int_t iSec, TVector3 &TLCS2GCSshift, TRotation &TLCS2GCS)
{
   TLCS2GCSshift = R0_tlsc[iSec];
   TLCS2GCS      = T_tl2glo[iSec];
}

//______________________________________________________________________

void TpcSectorGeoAlignmentVAK::GCS2TLCStransformation(const Int_t iSec, TVector3 &GCS2TLCSshift, TRotation &GCS2TLCS)
{
   GCS2TLCSshift = G0_tlsc[iSec];
   GCS2TLCS      = T_glo2tl[iSec];
}

//__________________________________________________________________

TVector2 TpcSectorGeoAlignmentVAK::LocalPadPosition(Int_t irow, Int_t ipad)
{
   /// Return local pad position for row and pad

   Double_t x, y;
   if (irow < fNrows[0]) {
      y = fPadH[0] * (irow + 0.5);
   } else {
      y = flYsec[1] - flYsec[0] + fPadH[1] * (irow - fNrows[0] + 0.5);
   }
   x = fPadW[1] * (ipad - fNPadsInRow[irow] / 2 + 0.5);
   return TVector2(x, y);
}

//__________________________________________________________________

void TpcSectorGeoAlignmentVAK::MinMaxPadPosition(Int_t irow, Int_t ipad, Double_t &xmin, Double_t &ymin, Double_t &xmax,
                                                 Double_t &ymax)
{
   /// Return local pad position for row and pad

   Double_t x, y;
   if (irow < fNrows[0]) {
      ymin = fPadH[0] * irow;
      ymax = ymin + fPadH[0];
   } else {
      ymin = flYsec[1] - flYsec[0] + fPadH[1] * (irow - fNrows[0]);
      ymax = ymin + fPadH[1];
   }
   xmin = fPadW[1] * (ipad - fNPadsInRow[irow] / 2);
   xmax = xmin + fPadW[1];
}

//__________________________________________________________________

TVector2 TpcSectorGeoAlignmentVAK::LocalPadPosition(Int_t padID)
{
   /// Return local pad position for padID

   Int_t    row = PadRow(padID);
   Int_t    pad = Pad(padID);
   Double_t x = 0.0, y = 0.0;
   if (row < fNrows[0]) {
      y = fPadH[0] * (row + 0.5);
   } else {
      y = flYsec[1] - flYsec[0] + fPadH[1] * (row - fNrows[0] + 0.5);
   }
   x = fPadW[1] * (pad - fNPadsInRow[row] / 2 + 0.5);
   return TVector2(x, y);
}

//__________________________________________________________________________

Int_t TpcSectorGeoAlignmentVAK::Global2Local(const Double_t *xyzGlob, Double_t *xyzLoc, Double_t &dist)
{
   /// Transform global coordinates to local (sector), calculate the distanse
   /// dist from the global point to the pad plane along the electrical field
   /// & return PadID for the projection

   /*
     Sector numbering for the camera is counterclockwise (looking at the end
   face outside the camera). For a camera with Z>0, the numbering starts at 0
   from the X axis (TPC coordinate system). For a camera with Z<0, the
   numbering starts at 12 from the negative x axis (TPC coordinate system).
     The coordinate axes of the local coordinate system of the sector form the
   right three. The Y axis is directed from the center of the chamber cylinder
   along it radius. The Z-axis of the sector for the positive part of the
   camera coincides with the direction of the z-axis of the global system,
   in the negative half they are opposite. From the point of view outside the
   camera on its end face, the direction of the local x axis is such that its
   rotation towards the local Y axis occurs counterclockwise for both cameras.
   The local coordinate center is in the center of the bottjv sector base
   (xGlob=yGlob=0) */

   /*The alignment perturbations are small, so the belonging of the projection
   point to the sector is determined for the canonical (ideal) position of the
   sector within the global coordinate system */

   Double_t xglob = xyzGlob[0];
   if (xyzGlob[2] < 0) xglob = -xglob;
   Double_t phGlob = TMath::ATan2(xyzGlob[1], xglob);
   Int_t    iSec;
   if (phGlob > 0) {
      iSec = (phGlob - fPhi0) / fDphi;
   } else {
      if (phGlob > fPhi0) {
         iSec = (phGlob - fPhi0) / fDphi;
      } else {
         iSec = (TMath::TwoPi() + phGlob - fPhi0) / fDphi;
      }
   }
   cout << "  phGlob; " << phGlob * TMath::RadToDeg() << endl;
   cout << "  fPhi0; " << fPhi0 * TMath::RadToDeg() << endl;
   if (xyzGlob[2] < 0) iSec += 12;
   Int_t result = iSec;
   // distance to the pads plane
   TVector3 vGlob(xyzGlob[0], xyzGlob[1], xyzGlob[2]);
   TVector3 vLoc = G0_gsc[iSec] + R_glo2loc[iSec] * vGlob;

   cout << "  iSec; " << iSec << endl;
   cout << vLoc.X() << "  " << vLoc.Y() << " " << vLoc.Z() << endl;
   cout << "  " << endl;
   cout << G0_gsc[iSec].X() << "  " << G0_gsc[iSec].Y() << " " << G0_gsc[iSec].Z() << endl;
   cout << "  " << endl;
   cout << R_glo2loc[iSec].XX() << "  " << R_glo2loc[iSec].XY() << " " << R_glo2loc[iSec].XZ() << endl;
   cout << R_glo2loc[iSec].YX() << "  " << R_glo2loc[iSec].YY() << " " << R_glo2loc[iSec].YZ() << endl;
   cout << R_glo2loc[iSec].ZX() << "  " << R_glo2loc[iSec].ZY() << " " << R_glo2loc[iSec].ZZ() << endl;
   cout << "  " << endl;
   cout << vGlob.X() << "  " << vGlob.Y() << " " << vGlob.Z() << endl;
   cout << "  " << endl;
   dist = vLoc.Z();
   // check the projection on pad plane is inside a pad
   if (vLoc.Y() < flYsec[0] || vLoc.Y() > flYsec[2]) { // outside the sector in Y
      result = -2;
   } else {
      Int_t irow = NofRow(vLoc.Y());      // carrent row number
      if (NofPad(irow, vLoc.X()) == -1) { // outside sector in X
         result = -1;
      }
   }
   xyzLoc[0] = vLoc.X();
   xyzLoc[1] = vLoc.Y();
   xyzLoc[2] = vLoc.Z();
   return result;
}

//__________________________________________________________________________

Int_t TpcSectorGeoAlignmentVAK::Global2Local(const TVector3 xyzGlob, TVector3 &xyzLoc, Double_t &dist)
{
   /// Transform global coordinates to local (sector)

   Double_t xyz1[3], xyz2[3], dist1;

   xyzGlob.GetXYZ(xyz1);
   Int_t iSec = Global2Local(xyz1, xyz2, dist1);
   xyzLoc.SetXYZ(xyz2[0], xyz2[1], xyz2[2]);
   dist = dist1;
   return iSec;
}

//__________________________________________________________________________

void TpcSectorGeoAlignmentVAK::Local2Global(Int_t iSec, const Double_t *xyzLoc, Double_t *xyzGlob)
{
   /// Transform local coordinates of sector iSec to global
   TVector3 Loc(xyzLoc[0], xyzLoc[1], xyzLoc[2]);
   TVector3 Glo = R0_lsc[iSec] + R_loc2glo[iSec] * Loc;
   xyzGlob[0]   = Glo.X();
   xyzGlob[1]   = Glo.Y();
   xyzGlob[2]   = Glo.Z();
   cout << "---l2g--- " << iSec << "------  " << endl;
   cout << R0_lsc[iSec].X() << "  " << R0_lsc[iSec].Y() << " " << R0_lsc[iSec].Z() << endl;
   cout << "  " << endl;
   cout << R_loc2glo[iSec].XX() << "  " << R_loc2glo[iSec].XY() << " " << R_loc2glo[iSec].XZ() << endl;
   cout << R_loc2glo[iSec].YX() << "  " << R_loc2glo[iSec].YY() << " " << R_loc2glo[iSec].YZ() << endl;
   cout << R_loc2glo[iSec].ZX() << "  " << R_loc2glo[iSec].ZY() << " " << R_loc2glo[iSec].ZZ() << endl;
   cout << "  " << endl;
}

//__________________________________________________________________________
void TpcSectorGeoAlignmentVAK::Local2Global(Int_t iSec, TVector3 xyzLoc, TVector3 &xyzGlob)
{
   /// Transform local coordinates of sector iSec to global
   Double_t xyz1[3], xyz2[3];
   xyzLoc.GetXYZ(xyz1);
   Local2Global(iSec, xyz1, xyz2);
   xyzGlob.SetXYZ(xyz2[0], xyz2[1], xyz2[2]);
}

//__________________________________________________________________________

Int_t TpcSectorGeoAlignmentVAK::NofRow(Double_t yloc)
{
   /// Compute row, pad and distance to the nearest row (+ or -)

   if (yloc < flYsec[0] || yloc > flYsec[2]) return -1;
   if (yloc <= flYsec[1]) {
      return Int_t(yloc / fPadH[0]);
   } else {
      return Int_t((yloc - flYsec[1]) / fPadH[1]) + fNrows[0];
   }
}

//__________________________________________________________________________

Bool_t TpcSectorGeoAlignmentVAK::IsItInside(Double_t xloc, Double_t yloc)
{
   /// check is the local point inside the sensitive sector area  (ROC)
   Int_t irow = NofRow(yloc);
   if (irow > 0) {
      Int_t    nPads = NofPadsInRow(irow);
      Double_t XPads = (nPads / 2) * fPadW[NofRoc(irow)];
      if (TMath::Abs(xloc) > XPads) {
         return false;
      } else {
         return true;
      }
   } else {
      return false;
   }
}

//__________________________________________________________________________

Bool_t TpcSectorGeoAlignmentVAK::IsItPad(Double_t xloc, Double_t yloc, Int_t &row, Int_t &pad, Double_t &x1,
                                         Double_t &x2, Double_t &y1, Double_t &y2)
{
   /// check is the local point inside the sensitive sector area  (ROC)
   /// return row & pad position numbers
   /// return minimun & maximun pad coordinates
   Int_t irow, ipad;
   irow = NofRow(yloc);
   if (irow >= 0) {
      Int_t    nPads = NofPadsInRow(irow);
      Double_t XPads = (nPads / 2) * fPadW[NofRoc(irow)];
      if (TMath::Abs(xloc) > XPads) {
         return false;
      } else {
         ipad = NofPad(irow, xloc);
         if (irow < fNrows[0]) {
            y1 = fPadH[0] * irow;
            y2 = y1 + fPadH[0];
         } else {
            y1 = flYsec[1] - flYsec[0] + fPadH[1] * (irow - fNrows[0]);
            y2 = y1 + fPadH[1];
         }
         x1  = fPadW[1] * (ipad - fNPadsInRow[irow] / 2);
         x2  = x1 + fPadW[1];
         row = irow;
         pad = ipad;
         return true;
      }
   } else {
      return false;
   }
}

//__________________________________________________________________________

void TpcSectorGeoAlignmentVAK::PadID(Double_t xloc, Double_t yloc, UInt_t &row, UInt_t &pad, Double_t &yNext)
{
   /// Compute row, pad and distance to the nearest row (+ or -)

   Double_t lowHeight = flYsec[1] - flYsec[0];
   if (yloc <= lowHeight) {
      row            = Int_t(yloc / fPadH[0]);
      pad            = Int_t(xloc / fPadW[0] + fNPadsInRow[row] / 2); // PadId is in the row
      Double_t delta = yloc - row * fPadH[0];
      if (delta > fPadH[0] / 2)
         yNext = fPadH[0] - delta;
      else
         yNext = -delta;
   } else {
      row            = Int_t((yloc - lowHeight) / fPadH[1]);
      pad            = Int_t(xloc / fPadW[1] + fNPadsInRow[fNrows[0] + row] / 2); // PadId is in the row
      Double_t delta = yloc - lowHeight - row * fPadH[0];
      if (delta > fPadH[1] / 2)
         yNext = fPadH[1] - delta;
      else
         yNext = -delta;
   }
}

//__________________________________________________________________________

Int_t TpcSectorGeoAlignmentVAK::NofPad(Int_t row, Double_t xloc)
{
   ///< get number of pad

   Double_t padW   = fPadW[NofRoc(row)];
   Int_t    nhpads = fNPadsInRow[row] / 2;
   if (TMath::Abs(xloc) >= padW * nhpads) return -1; // outside sector in X
   // pad number in the current row
   return Int_t((xloc + padW * nhpads) / padW);
}

//__________________________________________________________________________
// *************** obsolate **********************
Int_t TpcSectorGeoAlignmentVAK::Global2Local(const Double_t *xyzGlob, Double_t *xyzLoc,
                                             Int_t iSec) ///< transform global coordinates
                                                         /// to local (sector) - returns padID
{
   Double_t d;
   return Global2Local(xyzGlob, xyzLoc, d);
}

Int_t TpcSectorGeoAlignmentVAK::Global2Local(const TVector3 &xyzGlob, TVector3 &xyzLoc,
                                             Int_t iSec) ///< transform global coordinates to
                                                         ///  local (sector) - returns padID
{
   Double_t d;
   return Global2Local(xyzGlob, xyzLoc, d);
}

void TpcSectorGeoAlignmentVAK::PadID(Double_t xloc, Double_t yloc, UInt_t &row, UInt_t &pad, Float_t &yNext)
{
   /// Compute row, pad and distance to the nearest row (+ or -)

   Double_t lowHeight = flYsec[1] - flYsec[0];
   if (yloc <= lowHeight) {
      row            = Int_t(yloc / fPadH[0]);
      pad            = Int_t(xloc / fPadW[0] + fNPadsInRow[row] / 2); // PadId is in the row
      Double_t delta = yloc - row * fPadH[0];
      if (delta > fPadH[0] / 2)
         yNext = fPadH[0] - delta;
      else
         yNext = -delta + 0.000002;
   } else {
      row            = Int_t((yloc - lowHeight) / fPadH[1]);
      pad            = Int_t(xloc / fPadW[1] + fNPadsInRow[fNrows[0] + row] / 2); // PadId is in the row
      Double_t delta = yloc - lowHeight - row * fPadH[0];
      if (delta > fPadH[1] / 2)
         yNext = fPadH[1] - delta;
      else
         yNext = -delta + 0.000002;
   }
}
Int_t TpcSectorGeoAlignmentVAK::Global2Local(TVector3 xyzGlob, TVector3 &xyzLoc)
{
   Double_t d;
   return Global2Local(xyzGlob, xyzLoc, d);
}

// *************** obsolate **********************
