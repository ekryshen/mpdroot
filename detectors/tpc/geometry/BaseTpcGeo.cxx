//-----------------------------------------------------------
// Description:
//      BaseTpcGeo base class source file
//
// Environment:
//      Software for the MPD Detector at NICA.
//
// Authors:
//      Alexander Bychkov, Slavomir Hnatic
//      JINR, October, 2022
//-----------------------------------------------------------

#include "BaseTpcGeo.h"

//__________________________________________________________________________

BaseTpcGeo::BaseTpcGeo() {}

//__________________________________________________________________________

BaseTpcGeo::~BaseTpcGeo() {}

//__________________________________________________________________________

int BaseTpcGeo::SectorNumberFromGlobal(const TVector3 &globalXYZ)
{
   double phiGlobalNormalized = TMath::ATan2(globalXYZ.X(), globalXYZ.Y()) / SECTOR_PHI_RAD;
   // convention: 105 deg is start of Sector 0, 75 deg is start of Sector 1, ..
   const double PHI_SHIFT_NORMALIZED = 3.5;
   int          iSector;
   if (phiGlobalNormalized > PHI_SHIFT_NORMALIZED)
      iSector = static_cast<int>(PHI_SHIFT_NORMALIZED - phiGlobalNormalized + SECTOR_COUNT_HALF);
   else
      iSector = static_cast<int>(PHI_SHIFT_NORMALIZED - phiGlobalNormalized);

   if (globalXYZ.Z() < 0.0) iSector += SECTOR_COUNT_HALF;

   return iSector;
}

//__________________________________________________________________________

TVector2 BaseTpcGeo::PadRow2Local(double pad, double row)
{
   double x, y;
   int    rowInt = static_cast<int>(row);
   if (row > ROW_COUNT[inner]) {
      y = YPADAREA_LENGTH[inner] + (row - ROW_COUNT[inner]) * PAD_HEIGHT[outer];
      x = (pad - PAD_COUNT[rowInt]) * PAD_WIDTH[outer];
   } else {
      y = row * PAD_HEIGHT[inner];
      x = (pad - PAD_COUNT[rowInt]) * PAD_WIDTH[inner];
   }

   return {x, y};
}

//__________________________________________________________________________

TVector2 BaseTpcGeo::PadRowCenter2Local(int padNumber, int rowNumber)
{
   return PadRow2Local(padNumber + 0.5, rowNumber + 0.5);
}

//__________________________________________________________________________

std::pair<double, double> BaseTpcGeo::Local2PadRow(const TVector3 &localXYZ)
{
   double                          pad, row;
   const std::pair<double, double> ROW_OUTOFRANGE(0, -1.);
   const std::pair<double, double> PAD_OUTOFRANGE(1e4, 0);

   PadArea currentPadArea;
   if (localXYZ.Y() >= YPADAREA_LOCAL[lowerEdge] && localXYZ.Y() <= YPADAREA_LOCAL[midBoundary]) {
      currentPadArea = inner;
      row            = localXYZ.Y() / PAD_HEIGHT[inner];
   } else if (localXYZ.Y() > YPADAREA_LOCAL[midBoundary] && localXYZ.Y() < YPADAREA_LOCAL[upperEdge]) {
      currentPadArea = outer;
      row            = ROW_COUNT[inner] + (localXYZ.Y() - YPADAREA_LENGTH[inner]) / PAD_HEIGHT[outer];
   } else
      return ROW_OUTOFRANGE;

   pad = localXYZ.X() / PAD_WIDTH[currentPadArea] + PAD_COUNT[static_cast<int>(row)];

   if (pad < 0 || pad >= 2 * PAD_COUNT[static_cast<int>(row)]) return PAD_OUTOFRANGE;

   return std::make_pair(pad, row);
}

//__________________________________________________________________________

TVector3 BaseTpcGeo::Global2Local(const TVector3 &globalXYZ, int &iSector)
{
   iSector = SectorNumberFromGlobal(globalXYZ);

   TVector3 localXYZ(globalXYZ);

   localXYZ.RotateZ(SECTOR_PHI_RAD * iSector);

   if (localXYZ.Z() > 0) localXYZ.RotateY(TMath::Pi());

   localXYZ.SetY(localXYZ.Y() - YPADAREA_LOWEREDGE);
   localXYZ.SetZ(localXYZ.Z() + DRIFT_LENGTH);

   return localXYZ;
}

//__________________________________________________________________________

std::pair<TVector3, int> BaseTpcGeo::Global2Local(const TVector3 &globalXYZ)
{
   int      iSector;
   TVector3 localXYZ = Global2Local(globalXYZ, iSector);

   return std::make_pair(localXYZ, iSector);
}

//__________________________________________________________________________

TVector3 BaseTpcGeo::Local2Global(const TVector3 &localXYZ, int iSector)
{

   TVector3 globalXYZ(localXYZ);

   globalXYZ.SetZ(localXYZ.Z() - DRIFT_LENGTH);
   globalXYZ.SetY(localXYZ.Y() + YPADAREA_LOWEREDGE);

   if (iSector < SECTOR_COUNT_HALF) globalXYZ.RotateY(TMath::Pi());

   globalXYZ.RotateZ(-SECTOR_PHI_RAD * iSector);

   return globalXYZ;
}

//__________________________________________________________________________

TVector3 BaseTpcGeo::Local2Global(const std::pair<TVector3, int> &localPosition)
{
   return Local2Global(localPosition.first, localPosition.second);
}
