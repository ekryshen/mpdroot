#ifndef EIGENVECTOR_HH
#define EIGENVECTOR_HH
#include <vector>
#include "EigenCoefficients.h"
#include "Constants.h"
#include "Types.h"

namespace TpcAlignmentTest {
class EigenCalculatorTest;
class RunnerTest;
} // namespace TpcAlignmentTest

namespace TpcAlignmentLaserRays {
class Track;
class Point;
class CorrectionMatrix;

/// <summary>
/// Build eigen vector X=P*t+Q
/// Q  - center of mass of track points
/// P - major axis of the ellipsoid of inertia
/// </summary>
class EigenCalculator {
   friend class TpcAlignmentTest::EigenCalculatorTest;
   friend class TpcAlignmentTest::RunnerTest;

public:
   EigenCalculator()  = default;
   ~EigenCalculator() = default;

   /// <summary>
   /// Build Eigen Coefficients
   /// </summary>
   /// <param name="track">Input track</param>
   /// <param name="vCorrectionMatrix">A-matrix</param>
   /// <returns></returns>
   EigenCoefficients CalculateEigenValues(Track const &track, CorrectionMatrix const &vCorrectionMatrix);

   /// <summary>
   /// center of mass of track points
   /// </summary>
   /// <returns>Q-coefficient</returns>
   // Point const& Q() const;

   /// <summary>
   /// major axis of the ellipsoid of inertia
   /// </summary>
   /// <returns>P-coefficient</returns>
   // Point const& P() const;
protected:
   /// <summary>
   /// Calculate P-coefficient of track
   /// </summary>
   /// <param name="frontPoint"></param>
   /// <param name="backPoint"></param>
   /// <param name="vCorrectionMatrix">Correction matrix</param>
   /// <returns>P-coefficient</returns>
   Point CalculateP(Point const &frontPoint, Point const &backPoint, const CorrectionMatrix &vCorrectionMatrix) const;

   /// <summary>
   /// Calculate Q-coefficient
   /// </summary>
   /// <param name="track">Input track</param>
   /// <param name="vCorrectionMatrix">Correction matrix</param>

   Point CalculateQ(const Track &track, const CorrectionMatrix &vCorrectionMatrix);

   /// <summary>
   /// Accumulate coefficient to calculate Q
   /// </summary>
   /// <param name="TrackSize">number points on track</param>
   /// <param name="currentPoint">current point</param>
   /// <param name="vCorrectionMatrix">Correction matrix</param>
   static void AccumulateQ(Point &Q, double TrackSize, Point const &currentPoint,
                           const CorrectionMatrix &vCorrectionMatrix);

   /// <summary>
   /// Normilize coefficient P
   /// </summary>
   static void NormilizeP(Point &P);
};

} // namespace TpcAlignmentLaserRays

#endif
