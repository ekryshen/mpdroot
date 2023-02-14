#include "CorrectionMatrix.h"
#include "Linalg.h"
#include <cerrno>
#include <cstdio>
#include <iostream>
#include <cmath>
#include "FileHelper.h"
#include "Debug.h"
#include "Point.h"

namespace TpcAlignmentLaserRays {

CorrectionMatrix::CorrectionMatrix(/*SectorMatrix& a,*/ Debug const &inMode) : vDebug(inMode)
{
   vDebug.Print("testA start");
}

CorrectionMatrix::~CorrectionMatrix()
{
   vDebug.Print("testA finished");
}

SectorMatrix const &CorrectionMatrix::Get() const
{
   return vCorrectionMatrix;
}

SectorMatrix &CorrectionMatrix::Get()
{
   return vCorrectionMatrix;
}

void CorrectionMatrix::SetMatrixWithDefaultValues()
{
   vDebug.Print("SetCorrectionMatrix2Default\n");
   for (int k = 0; k < SECTORS; k++) {
      vCorrectionMatrix.at(k).resize(L_);
      for (size_t i = 0; i < 3u; i++) {    // i - point
         for (size_t j = 0; j < M_; j++) { // j - point's coordinate
            vCorrectionMatrix[k][i + 3u * j] = 0.;
         }
         vCorrectionMatrix[k][i + 3 * i] = 1.;
      }
   }
}

void CorrectionMatrix::LoadMatrix4File(std::string filename)
{
   if (filename.empty()) {
      throw std::runtime_error("Correction Matrix file path is not set.");
   }
   FileHelper::LoadCorrectionMatrix(filename, vCorrectionMatrix);
}

void CorrectionMatrix::SaveMatrix2File(std::string filename, int precision)
{
   if (filename.empty()) {
      throw std::runtime_error("Correction Matrix file path is not set.");
   }
   FileHelper::SaveCorrectionMatrix(filename, vCorrectionMatrix, precision);
}

double CorrectionMatrix::CalculateCoordinatesNewValue(Point const &point, unsigned coordinate) const
{
   size_t i{coordinate};

   return vCorrectionMatrix[point.Sector()][i + static_cast<size_t>(3 * 0)] * point.X    // i+3*0
          + vCorrectionMatrix[point.Sector()][i + static_cast<size_t>(3 * 1)] * point.Y  // i+3*1
          + vCorrectionMatrix[point.Sector()][i + static_cast<size_t>(3 * 2)] * point.Z; // i+3*2
}

Point CorrectionMatrix::CorrectPoint(Point const &inputPoint) const
{
   Point result{CalculateCoordinatesNewValue(inputPoint, 0u), CalculateCoordinatesNewValue(inputPoint, 1u),
                CalculateCoordinatesNewValue(inputPoint, 2u)};
   result.SetSector();
   return result;
}

} // namespace TpcAlignmentLaserRays
