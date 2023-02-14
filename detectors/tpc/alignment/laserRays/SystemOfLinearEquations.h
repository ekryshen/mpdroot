#ifndef SYSTEMOFLINEAREQUATIONS_HH
#define SYSTEMOFLINEAREQUATIONS_HH
#include <vector>
#include <array>
#include "Constants.h"
#include "Types.h"

namespace TpcAlignmentTest {
class SystemOfLinearEquationsTest;
class RunnerTest;
} // namespace TpcAlignmentTest

namespace TpcAlignmentLaserRays {
class Rays;
class Track;
class Point;
class CorrectionMatrix;
struct EigenCoefficients;

/// <summary>
/// Operate on System of Linear Equations
/// </summary>
class SystemOfLinearEquations {
   friend class TpcAlignmentTest::SystemOfLinearEquationsTest;
   friend class TpcAlignmentTest::RunnerTest;

private:
   /// <summary>
   /// Left side of equation
   /// </summary>
   SectorMatrix vM;

   /// <summary>
   /// Right side of equation
   /// </summary>
   SectorMatrix vR;

public:
   SystemOfLinearEquations();
   ~SystemOfLinearEquations() = default;
   /// <summary>
   /// Build and Solve System of Linear Equations
   /// </summary>
   /// <param name="tracks">input Rays</param>
   /// <param name="A">A-correction matrice</param>
   void BuildAndSolve(Rays &tracks, CorrectionMatrix &A);

   /// <summary>
   /// Get M-coefficients of the given sector
   /// </summary>
   /// <param name="sector">Sector index</param>
   /// <returns>M-coefficients</returns>
   std::vector<double> const &M(int sector) const { return vM[sector]; }

   /// <summary>
   /// Get R-coefficients of the given sector
   /// </summary>
   /// <param name="sector">Sector index</param>
   /// <returns>R-coefficients</returns>
   std::vector<double> const &R(int sector) const { return vR[sector]; }

   SectorMatrix const &GetM() const { return vM; }

   SectorMatrix const &GetR() const { return vR; }

private:
   /// <summary>
   /// Utility function to AccumulateMR
   /// </summary>
   /// <param name="track"></param>
   /// <param name="eigenCoefficients"></param>
   void AccumulateMR(const Track &track, const EigenCoefficients &eigenCoefficients);
   /// <summary>
   /// Calculate G-coeff
   /// </summary>
   static std::vector<double> CalculateG(Point const &P, Point const &point, double cf = 1.);
   /// <summary>
   /// Calculate H-coeff
   /// </summary>
   static std::vector<double> CalculateH(Point const &P, Point const &Q, Point const &point, double cf = 1.);
   /// <summary>
   /// Calculate MR-coeff for a sector
   /// </summary>
   void CalculateMR(int sector, std::vector<double> const &G, std::vector<double> const &H);
   /// <summary>
   /// Solve AMR equations - sector by sector
   /// </summary>
   void SolveAmrEqueations(SectorMatrix &A);

   /// <summary>
   /// Build M&R-coeff
   /// </summary>
   void BuildMR(Rays const &tracks, CorrectionMatrix const &A);
};
} // namespace TpcAlignmentLaserRays
#endif
