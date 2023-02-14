#include "EigenCalculator.h"
#include "Track.h"
#include "Point.h"
#include "EigenCoefficients.h"
#include <stdexcept>
#include "Constants.h"
#include "CorrectionMatrix.h"

namespace TpcAlignmentLaserRays {

EigenCoefficients EigenCalculator::CalculateEigenValues(Track const &track, CorrectionMatrix const &vCorrectionMatrix)
{
   if (track.Get().empty()) {
      return {{}, {}};
   }

   Point Q{CalculateQ(track, vCorrectionMatrix)};
   Point P{CalculateP(track.Get().front(), track.Get().back(), vCorrectionMatrix)};

   return {P, Q};
}

Point EigenCalculator::CalculateQ(const Track &track, const CorrectionMatrix &vCorrectionMatrix)
{
   Point        Q{0, 0, 0};
   double const vTrackSize{static_cast<double>(track.Get().size())};
   for (size_t i = 0; i < track.Get().size(); ++i) {
      AccumulateQ(Q, vTrackSize, track.Get().at(i), vCorrectionMatrix);
   }
   return Q;
}

Point EigenCalculator::CalculateP(Point const &frontPoint, Point const &backPoint,
                                  const CorrectionMatrix &vCorrectionMatrix) const
{
   Point back  = vCorrectionMatrix.CorrectPoint(backPoint);
   Point front = vCorrectionMatrix.CorrectPoint(frontPoint);
   Point P     = back - front;
   NormilizeP(P);

   return P;
}

void EigenCalculator::NormilizeP(Point &P)
{
   double const r{sqrt(P.X * P.X + P.Y * P.Y + P.Z * P.Z)};
   double const tolerance{1.E-32};
   if (fabs(r) < tolerance) {
      throw std::runtime_error("R is quite close to Zero");
   }
   P /= r;
}

void EigenCalculator::AccumulateQ(Point &Q, double TrackSize, const Point &currentPoint,
                                  const CorrectionMatrix &vCorrectionMatrix)
{
   Point correctedPoint{vCorrectionMatrix.CorrectPoint(currentPoint)};
   correctedPoint /= TrackSize;
   Q += correctedPoint;
}
} // namespace TpcAlignmentLaserRays
