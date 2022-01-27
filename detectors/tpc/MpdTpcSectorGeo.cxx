/// \class MpdTpcSectorGeo
///
/// MPD TPC sector geometry description
/// \author Alexander Zinchenko (LHEP, JINR, Dubna)

#include "MpdTpcSectorGeo.h"
#include "TpcGas.h"
//#include "TpcGeoPar.h"
//#include "FairGeoNode.h"
//#include "FairRunAna.h"
//#include "FairRuntimeDb.h"

#include <TGeoManager.h>
#include <TGeoPgon.h>
#include <TGeoTube.h>
#include <TMath.h>
#include <TSystem.h>
#include <Riostream.h>

using std::cout;
using std::endl;

MpdTpcSectorGeo *MpdTpcSectorGeo::fgTpcSec = 0x0;

//__________________________________________________________________________
MpdTpcSectorGeo *MpdTpcSectorGeo::Instance()
{
   /// Get pointer to the TPC sector singleton object
   if (!fgTpcSec) {
      fgTpcSec = new MpdTpcSectorGeo;

      // Automatic deletion
      std::atexit(DestroyInstance);
   }
   return fgTpcSec;
}

//__________________________________________________________________________
void MpdTpcSectorGeo::Init()
{
   /// Initialization of sector parameters

   // FairRuntimeDb* rtdb = FairRun::Instance()->GetRuntimeDb();
   // TpcGeoPar *geoPar = (TpcGeoPar*) rtdb->getContainer("TpcGeoPar");
   TString     volName = "tpc01sv";
   TGeoVolume *sv      = gGeoManager->GetVolume(volName);
   if (!sv) {
      Warning("Init", "No sensitive volume found !!!");
      return;
   }
   TGeoShape *shape = sv->GetShape();
   // TObjArray* sensNodes = geoPar->GetGeoSensitiveNodes();
   // cout << sensNodes->GetEntriesFast() << " " << geoPar->GetGeoPassiveNodes()->GetEntriesFast() << endl;
   // FairGeoNode* sensVol = (FairGeoNode*) (sensNodes->FindObject(volName));
   // FairGeoNode* sensVol = (FairGeoNode*) (sensNodes->At(0));
   // TArrayD* params = sensVol->getParameters();
   fDphi = TMath::TwoPi() / fgkNsect;
   fPhi0 = -fDphi / 2.;

   cout << "\n !!!!! ***** ***** ***** !!!!!" << endl;
   // if (sensVol->getShape() == "TUBE") {
   if (TString(shape->ClassName()) == "TGeoTube") {
      // Old sensitive volume
      fNrows[0] = fNrows[1] = 50;
      // fNrows = 53;
      // Double_t rMin = params->At(0);
      // Double_t rMax = params->At(1);
      Double_t rMin = ((TGeoTube *)shape)->GetRmin();
      Double_t rMax = ((TGeoTube *)shape)->GetRmax();
      cout << " ***** TPC sensitive volume: " << sv->GetName() << " " << shape->ClassName() << " " << rMin << " "
           << rMax << " " << ((TGeoTube *)shape)->GetDZ() << endl;
      // fYsec[0] = rMin;
      // fYsec[1] = rMax * TMath::Cos(fPhi0);
      fYsec[0] = 36.;            // for consistency with earlier results
      fYsec[1] = fYsec[2] = 96.; // for consistency with earlier results
      fPadH[0] = fPadH[1] = (fYsec[1] - fYsec[0]) / fNrows[0];
      // fYsec[0] = 40.4;
      // fYsec[1] = fYsec[0] + 53.0;
      //} else if (sensVol->getShape() == "PGON") {
   } else if (TString(shape->ClassName()) == "TGeoPgon") {
      // New sensitive volume
      // fYsec[0] = params->At(5);
      // fYsec[1] = params->At(6);
      fYsec[0] = ((TGeoPgon *)shape)->Rmin(0) + 0.4;
      // fYsec[2] = ((TGeoPgon*)shape)->Rmax(0);
      fPadH[0]    = 1.2; // approx. pad height in inner ROC region
      fPadH[1]    = 1.8; // 1.2; // approx. pad height in outer ROC region
      fNrows[0]   = 27;  // 33;
      fNrows[1]   = 26;  // 33;
      Double_t dy = fPadH[0] * fNrows[0] + fPadH[1] * fNrows[1];
      fYsec[2]    = fYsec[0] + dy;
      // Double_t scale = (fYsec[2] - fYsec[0]) / dy;
      // fPadH[0] *= scale;
      // fPadH[1] *= scale;
      fYsec[1] = fYsec[0] + fPadH[0] * fNrows[0];
      // fYsec[1] = fYsec[0] + 60;
      // fNrows = Int_t ((fYsec[1] - fYsec[0] + 0.1) / 1.2); // 66 lays
      // fNrows = Int_t ((fYsec[1] - fYsec[0] + 0.1) / 1.5); // 53 lays
      // fNrows = Int_t ((fYsec[1] - fYsec[0] + 0.1) / 1.0);
      // Starting phi
      fPhi0           = -((TGeoPgon *)shape)->Phi1() * TMath::DegToRad();
      Double_t loc[3] = {100, 0, 10}, glob[3] = {0};
      gGeoManager->FindNode(loc[0], loc[1], loc[2]);
      gGeoManager->LocalToMaster(loc, glob);
      fPhi0 += TMath::ATan2(glob[1], glob[0]); // due to rotation
      // Extract sensitive volume Z-coordinates
      TGeoVolume *membr = gGeoManager->GetVolume("tpc01mb"); // membrane
      fZmin             = ((TGeoTube *)membr->GetShape())->GetDZ();
      // fZmax = fZmin + 2 * ((TGeoTube*)shape)->GetDZ();
      fZmax = 170.;
      cout << " ***** TPC sensitive volume: " << sv->GetName() << ", shape: " << shape->ClassName()
           << ", inner radius: " << fYsec[0] << ", outer radius: " << fYsec[2] << ", Zmin, Zmax: " << fZmin << ", "
           << fZmax << endl;
   } else
      Fatal("MpdTpcSectorGeo::Init()", " !!! Unknown sensitive volume shape !!! ");

   // fPadH = (fYsec[1] - fYsec[0]) / fNrows;
   cout << " ***** TPC sector params - inner radius: " << fYsec[0] << ", outer radius: " << fYsec[2]
        << ", boundary: " << fYsec[1] << endl;
   cout << " Number of sectors: " << fgkNsect << ", phi0: " << fPhi0 * TMath::RadToDeg()
        << ", numbers of padrows: " << fNrows[0] << " " << fNrows[1] << ", pad heights: " << fPadH[0] << " " << fPadH[1]
        << endl;
   cout << " !!!!! ***** ***** ***** !!!!!" << endl;

   fPadW[0] = fPadW[1] = 0.5; // pad widths
   // Numbers of pads in rows
   // Double_t tan = TMath::Tan(fDphi/2), dead = 1.35; // dead area on one side
   // Double_t tan = TMath::Tan(fDphi/2), dead = 0.15; // dead area on one side
   fNPadsInRows = new Int_t[NofRows()];
   // for (Int_t j = 0; j < fNrows[0]; ++j)
   //  fNPadsInRows[j] = Int_t ((tan * (fYsec[0] + j * fPadH[0]) - dead) / fPadW[0]);
   // for (Int_t j = 0; j < fNrows[1]; ++j)
   //  fNPadsInRows[j+fNrows[0]] = Int_t ((tan * (fYsec[1] + j * fPadH[1]) - dead) / fPadW[1]);
   Int_t nPiR[] = {20, 21, 21, 22, 23, 23, 24, 24, 25, 26, 26, 27, 28, 28, 29, 30, 30, 31,
                   32, 32, 33, 33, 34, 35, 35, 36, 37, 38, 39, 40, 41, 41, 42, 43, 44, 45,
                   46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62};
   for (Int_t j = 0; j < NofRows(); ++j) fNPadsInRows[j] = nPiR[j];

   // Gas parameters
   std::string tpcGasFile = gSystem->Getenv("VMCWORKDIR");
   tpcGasFile += "/geometry/Ar-90_CH4-10.asc";
   fGas     = new TpcGas(tpcGasFile, 130);
   fTimeMax = fZmax / fGas->VDrift();

   fNTimeBins = 512; // max number of time bins
   // fTimeBin = fTimeMax / fNTimeBins;
   fTimeBin = 100; // 100 ns
   // fZ2TimeBin = fNTimeBins / fZmax;
   // fTimeBinMax = Int_t (fTimeMax / fTimeBin);
   fTimeBinMax = fTimeMax / fTimeBin;
   fZ2TimeBin  = fTimeBinMax / fZmax;
}

//__________________________________________________________________________
Int_t MpdTpcSectorGeo::Global2Local(const Double_t *xyzGlob, Double_t *xyzLoc, Int_t iSec0)
{
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
   if (x < fYsec[0] + safety || x > fYsec[2] - safety) return -1; // outside sector in Y
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
Int_t MpdTpcSectorGeo::Global2Local(const TVector3 &xyzGlob, TVector3 &xyzLoc, Int_t iSec0)
{
   /// Transform global coordinates to local (sector)

   Double_t xyz1[3], xyz2[3];

   xyzGlob.GetXYZ(xyz1);
   Int_t iSec = Global2Local(xyz1, xyz2, iSec0);
   xyzLoc.SetXYZ(xyz2[0], xyz2[1], xyz2[2]);
   return iSec;
}

//__________________________________________________________________________
void MpdTpcSectorGeo::Local2Global(Int_t iSec, const Double_t *xyzLoc, Double_t *xyzGlob)
{
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
void MpdTpcSectorGeo::Local2Global(Int_t iSec, const TVector3 &xyzLoc, TVector3 &xyzGlob)
{
   /// Transform local coordinates of sector iSec to global

   Double_t xyz1[3], xyz2[3];

   xyzLoc.GetXYZ(xyz1);
   Local2Global(iSec, xyz1, xyz2);
   xyzGlob.SetXYZ(xyz2[0], xyz2[1], xyz2[2]);
}

//__________________________________________________________________________
TVector2 MpdTpcSectorGeo::LocalPadPosition(Int_t padID)
{
   /// Return local pad position for padID

   Int_t    row = PadRow(padID);
   Double_t x = 0.0, y = 0.0;
   if (row < fNrows[0])
      y = fPadH[0] * (row + 0.5);
   else
      y = fYsec[1] - fYsec[0] + fPadH[1] * (row - fNrows[0] + 0.5);
   return TVector2(x, y);
}

//__________________________________________________________________________

void MpdTpcSectorGeo::PadID(Float_t xloc, Float_t yloc, UInt_t &row, UInt_t &pad, Float_t &yNext)
{
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

//__________________________________________________________________________
