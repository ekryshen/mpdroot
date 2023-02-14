#include "SystemOfLinearEquations.h"
#include "Track.h"
#include "Rays.h"
#include "Point.h"
#include <vector>
#include "EigenCoefficients.h"
#include "Linalg.h"
#include "CorrectionMatrix.h"
#include "EigenCalculator.h"

using std::vector;

namespace TpcAlignmentLaserRays {

SystemOfLinearEquations::SystemOfLinearEquations()
{
   for (size_t i = 0; i < SECTORS; i++) {
      vM[i].resize(LL, 0.0);
      vR[i].resize(L_, 0.0);
   }
}

void SystemOfLinearEquations::AccumulateMR(const Track &track, const EigenCoefficients &eigenCoefficients)
{
   for (size_t i = 0; i < track.Get().size(); ++i) {
      CalculateMR(track.Get().at(i).Sector(), CalculateG(eigenCoefficients.P, track.Get().at(i)),
                  CalculateH(eigenCoefficients.P, eigenCoefficients.Q, track.Get().at(i)));
   }
}

void SystemOfLinearEquations::CalculateMR(int sector, vector<double> const &G, vector<double> const &H)
{
   for (size_t i = 0; i < L_; i++) {
      for (size_t j = 0; j < L_; j++) {
         // M[sector * 144 + i + j * 12] = M[sector * 144 + i + j * 12] + G[i + j * 12];
         vM[sector][i + j * L_] += G[i + j * L_];
      }
      vR[sector][i] -= H[i];
      // R[sector * 12 + i] = R[sector * 12 + i] - H[i];
   }
}

void SystemOfLinearEquations::BuildMR(Rays const &tracks, CorrectionMatrix const &A)
{
   for (Track const &track : tracks.Get()) {
      EigenCalculator   vEigenValues;
      EigenCoefficients coeff = vEigenValues.CalculateEigenValues(track, A);
      AccumulateMR(track, coeff);
   }
}

void SystemOfLinearEquations::BuildAndSolve(Rays &tracks, CorrectionMatrix &A)
{
   BuildMR(tracks, A);
   SolveAmrEqueations(A.Get());
}

void SystemOfLinearEquations::SolveAmrEqueations(SectorMatrix &A)
{
   for (size_t i = 0; i < SECTORS; i++) {
      Linalg::lsolve(vM[i], vR[i], A[i]);
   }
}

vector<double> SystemOfLinearEquations::CalculateG(Point const &P, Point const &point, double cf)
{
   {
      vector<double> G(144, 0.0);
      G[0]   = P.Y * P.Y * point.X * point.X + P.Z * P.Z * point.X * point.X;
      G[12]  = -P.X * point.X * point.X * P.Y;
      G[24]  = -P.X * point.X * point.X * P.Z;
      G[36]  = P.Y * P.Y * point.X * point.Y + P.Z * P.Z * point.X * point.Y;
      G[48]  = -P.X * point.Y * P.Y * point.X;
      G[60]  = -P.X * point.Y * P.Z * point.X;
      G[72]  = P.Y * P.Y * point.X * point.Z + P.Z * P.Z * point.X * point.Z;
      G[84]  = -P.X * point.Z * P.Y * point.X;
      G[96]  = -P.X * point.Z * P.Z * point.X;
      G[108] = P.Y * P.Y * point.X * cf + P.Z * P.Z * point.X * cf;
      G[120] = -P.X * cf * P.Y * point.X;
      G[132] = -P.X * cf * P.Z * point.X;
      G[1]   = -P.X * point.X * point.X * P.Y;
      G[13]  = P.X * P.X * point.X * point.X + P.Z * P.Z * point.X * point.X;
      G[25]  = -P.Y * point.X * point.X * P.Z;
      G[37]  = -P.X * point.Y * P.Y * point.X;
      G[49]  = P.X * P.X * point.X * point.Y + P.Z * P.Z * point.X * point.Y;
      G[61]  = -P.Y * point.Y * P.Z * point.X;
      G[73]  = -P.X * point.Z * P.Y * point.X;
      G[85]  = P.X * P.X * point.X * point.Z + P.Z * P.Z * point.X * point.Z;
      G[97]  = -P.Y * point.Z * P.Z * point.X;
      G[109] = -P.X * cf * P.Y * point.X;
      G[121] = P.X * P.X * point.X * cf + P.Z * P.Z * point.X * cf;
      G[133] = -P.Y * cf * P.Z * point.X;
      G[2]   = -P.X * point.X * point.X * P.Z;
      G[14]  = -P.Y * point.X * point.X * P.Z;
      G[26]  = P.X * P.X * point.X * point.X + P.Y * P.Y * point.X * point.X;
      G[38]  = -P.X * point.Y * P.Z * point.X;
      G[50]  = -P.Y * point.Y * P.Z * point.X;
      G[62]  = P.X * P.X * point.X * point.Y + P.Y * P.Y * point.X * point.Y;
      G[74]  = -P.X * point.Z * P.Z * point.X;
      G[86]  = -P.Y * point.Z * P.Z * point.X;
      G[98]  = P.X * P.X * point.X * point.Z + P.Y * P.Y * point.X * point.Z;
      G[110] = -P.X * cf * P.Z * point.X;
      G[122] = -P.Y * cf * P.Z * point.X;
      G[134] = P.X * P.X * point.X * cf + P.Y * P.Y * point.X * cf;
      G[3]   = P.Y * P.Y * point.X * point.Y + P.Z * P.Z * point.X * point.Y;
      G[15]  = -P.X * point.Y * P.Y * point.X;
      G[27]  = -P.X * point.Y * P.Z * point.X;
      G[39]  = P.Y * P.Y * point.Y * point.Y + P.Z * P.Z * point.Y * point.Y;
      G[51]  = -P.X * point.Y * point.Y * P.Y;
      G[63]  = -P.X * point.Y * point.Y * P.Z;
      G[75]  = P.Y * P.Y * point.Y * point.Z + P.Z * P.Z * point.Y * point.Z;
      G[87]  = -P.X * point.Z * P.Y * point.Y;
      G[99]  = -P.X * point.Z * P.Z * point.Y;
      G[111] = P.Y * P.Y * point.Y * cf + P.Z * P.Z * point.Y * cf;
      G[123] = -P.X * cf * P.Y * point.Y;
      G[135] = -P.X * cf * P.Z * point.Y;
      G[4]   = -P.X * point.Y * P.Y * point.X;
      G[16]  = P.X * P.X * point.X * point.Y + P.Z * P.Z * point.X * point.Y;
      G[28]  = -P.Y * point.Y * P.Z * point.X;
      G[40]  = -P.X * point.Y * point.Y * P.Y;
      G[52]  = P.X * P.X * point.Y * point.Y + P.Z * P.Z * point.Y * point.Y;
      G[64]  = -P.Y * point.Y * point.Y * P.Z;
      G[76]  = -P.X * point.Z * P.Y * point.Y;
      G[88]  = P.X * P.X * point.Y * point.Z + P.Z * P.Z * point.Y * point.Z;
      G[100] = -P.Y * point.Z * P.Z * point.Y;
      G[112] = -P.X * cf * P.Y * point.Y;
      G[124] = P.X * P.X * point.Y * cf + P.Z * P.Z * point.Y * cf;
      G[136] = -P.Y * cf * P.Z * point.Y;
      G[5]   = -P.X * point.Y * P.Z * point.X;
      G[17]  = -P.Y * point.Y * P.Z * point.X;
      G[29]  = P.X * P.X * point.X * point.Y + P.Y * P.Y * point.X * point.Y;
      G[41]  = -P.X * point.Y * point.Y * P.Z;
      G[53]  = -P.Y * point.Y * point.Y * P.Z;
      G[65]  = P.X * P.X * point.Y * point.Y + P.Y * P.Y * point.Y * point.Y;
      G[77]  = -P.X * point.Z * P.Z * point.Y;
      G[89]  = -P.Y * point.Z * P.Z * point.Y;
      G[101] = P.X * P.X * point.Y * point.Z + P.Y * P.Y * point.Y * point.Z;
      G[113] = -P.X * cf * P.Z * point.Y;
      G[125] = -P.Y * cf * P.Z * point.Y;
      G[137] = P.X * P.X * point.Y * cf + P.Y * P.Y * point.Y * cf;
      G[6]   = P.Y * P.Y * point.X * point.Z + P.Z * P.Z * point.X * point.Z;
      G[18]  = -P.X * point.Z * P.Y * point.X;
      G[30]  = -P.X * point.Z * P.Z * point.X;
      G[42]  = P.Y * P.Y * point.Y * point.Z + P.Z * P.Z * point.Y * point.Z;
      G[54]  = -P.X * point.Z * P.Y * point.Y;
      G[66]  = -P.X * point.Z * P.Z * point.Y;
      G[78]  = P.Y * P.Y * point.Z * point.Z + P.Z * P.Z * point.Z * point.Z;
      G[90]  = -P.X * point.Z * point.Z * P.Y;
      G[102] = -P.X * point.Z * point.Z * P.Z;
      G[114] = P.Y * P.Y * point.Z * cf + P.Z * P.Z * point.Z * cf;
      G[126] = -P.X * cf * P.Y * point.Z;
      G[138] = -P.X * cf * P.Z * point.Z;
      G[7]   = -P.X * point.Z * P.Y * point.X;
      G[19]  = P.X * P.X * point.X * point.Z + P.Z * P.Z * point.X * point.Z;
      G[31]  = -P.Y * point.Z * P.Z * point.X;
      G[43]  = -P.X * point.Z * P.Y * point.Y;
      G[55]  = P.X * P.X * point.Y * point.Z + P.Z * P.Z * point.Y * point.Z;
      G[67]  = -P.Y * point.Z * P.Z * point.Y;
      G[79]  = -P.X * point.Z * point.Z * P.Y;
      G[91]  = P.X * P.X * point.Z * point.Z + P.Z * P.Z * point.Z * point.Z;
      G[103] = -P.Y * point.Z * point.Z * P.Z;
      G[115] = -P.X * cf * P.Y * point.Z;
      G[127] = P.X * P.X * point.Z * cf + P.Z * P.Z * point.Z * cf;
      G[139] = -P.Y * cf * P.Z * point.Z;
      G[8]   = -P.X * point.Z * P.Z * point.X;
      G[20]  = -P.Y * point.Z * P.Z * point.X;
      G[32]  = P.X * P.X * point.X * point.Z + P.Y * P.Y * point.X * point.Z;
      G[44]  = -P.X * point.Z * P.Z * point.Y;
      G[56]  = -P.Y * point.Z * P.Z * point.Y;
      G[68]  = P.X * P.X * point.Y * point.Z + P.Y * P.Y * point.Y * point.Z;
      G[80]  = -P.X * point.Z * point.Z * P.Z;
      G[92]  = -P.Y * point.Z * point.Z * P.Z;
      G[104] = P.X * P.X * point.Z * point.Z + P.Y * P.Y * point.Z * point.Z;
      G[116] = -P.X * cf * P.Z * point.Z;
      G[128] = -P.Y * cf * P.Z * point.Z;
      G[140] = P.X * P.X * point.Z * cf + P.Y * P.Y * point.Z * cf;
      G[9]   = P.Y * P.Y * point.X * cf + P.Z * P.Z * point.X * cf;
      G[21]  = -P.X * cf * P.Y * point.X;
      G[33]  = -P.X * cf * P.Z * point.X;
      G[45]  = P.Y * P.Y * point.Y * cf + P.Z * P.Z * point.Y * cf;
      G[57]  = -P.X * cf * P.Y * point.Y;
      G[69]  = -P.X * cf * P.Z * point.Y;
      G[81]  = P.Y * P.Y * point.Z * cf + P.Z * P.Z * point.Z * cf;
      G[93]  = -P.X * cf * P.Y * point.Z;
      G[105] = -P.X * cf * P.Z * point.Z;
      G[117] = P.Y * P.Y * cf * cf + P.Z * P.Z * cf * cf;
      G[129] = -P.X * cf * cf * P.Y;
      G[141] = -P.X * cf * cf * P.Z;
      G[10]  = -P.X * cf * P.Y * point.X;
      G[22]  = P.X * P.X * point.X * cf + P.Z * P.Z * point.X * cf;
      G[34]  = -P.Y * cf * P.Z * point.X;
      G[46]  = -P.X * cf * P.Y * point.Y;
      G[58]  = P.X * P.X * point.Y * cf + P.Z * P.Z * point.Y * cf;
      G[70]  = -P.Y * cf * P.Z * point.Y;
      G[82]  = -P.X * cf * P.Y * point.Z;
      G[94]  = P.X * P.X * point.Z * cf + P.Z * P.Z * point.Z * cf;
      G[106] = -P.Y * cf * P.Z * point.Z;
      G[118] = -P.X * cf * cf * P.Y;
      G[130] = P.X * P.X * cf * cf + P.Z * P.Z * cf * cf;
      G[142] = -P.Y * cf * cf * P.Z;
      G[11]  = -P.X * cf * P.Z * point.X;
      G[23]  = -P.Y * cf * P.Z * point.X;
      G[35]  = P.X * P.X * point.X * cf + P.Y * P.Y * point.X * cf;
      G[47]  = -P.X * cf * P.Z * point.Y;
      G[59]  = -P.Y * cf * P.Z * point.Y;
      G[71]  = P.X * P.X * point.Y * cf + P.Y * P.Y * point.Y * cf;
      G[83]  = -P.X * cf * P.Z * point.Z;
      G[95]  = -P.Y * cf * P.Z * point.Z;
      G[107] = P.X * P.X * point.Z * cf + P.Y * P.Y * point.Z * cf;
      G[119] = -P.X * cf * cf * P.Z;
      G[131] = -P.Y * cf * cf * P.Z;
      G[143] = P.X * P.X * cf * cf + P.Y * P.Y * cf * cf;
      return G;
   }
}

vector<double> SystemOfLinearEquations::CalculateH(Point const &P, Point const &Q, Point const &point, double cf)
{
   vector<double> H(12, 0.0);
   H[0]  = -(-P.X * Q.Z + P.Z * Q.X) * P.Z * point.X + (P.X * Q.Y - P.Y * Q.X) * P.Y * point.X;
   H[1]  = (P.Y * Q.Z - P.Z * Q.Y) * P.Z * point.X - (P.X * Q.Y - P.Y * Q.X) * P.X * point.X;
   H[2]  = -(P.Y * Q.Z - P.Z * Q.Y) * P.Y * point.X + (-P.X * Q.Z + P.Z * Q.X) * P.X * point.X;
   H[3]  = -(-P.X * Q.Z + P.Z * Q.X) * P.Z * point.Y + (P.X * Q.Y - P.Y * Q.X) * P.Y * point.Y;
   H[4]  = (P.Y * Q.Z - P.Z * Q.Y) * P.Z * point.Y - (P.X * Q.Y - P.Y * Q.X) * P.X * point.Y;
   H[5]  = -(P.Y * Q.Z - P.Z * Q.Y) * P.Y * point.Y + (-P.X * Q.Z + P.Z * Q.X) * P.X * point.Y;
   H[6]  = -(-P.X * Q.Z + P.Z * Q.X) * P.Z * point.Z + (P.X * Q.Y - P.Y * Q.X) * P.Y * point.Z;
   H[7]  = (P.Y * Q.Z - P.Z * Q.Y) * P.Z * point.Z - (P.X * Q.Y - P.Y * Q.X) * P.X * point.Z;
   H[8]  = -(P.Y * Q.Z - P.Z * Q.Y) * P.Y * point.Z + (-P.X * Q.Z + P.Z * Q.X) * P.X * point.Z;
   H[9]  = -(-P.X * Q.Z + P.Z * Q.X) * P.Z * cf + (P.X * Q.Y - P.Y * Q.X) * P.Y * cf;
   H[10] = (P.Y * Q.Z - P.Z * Q.Y) * P.Z * cf - (P.X * Q.Y - P.Y * Q.X) * P.X * cf;
   H[11] = -(P.Y * Q.Z - P.Z * Q.Y) * P.Y * cf + (-P.X * Q.Z + P.Z * Q.X) * P.X * cf;
   return H;
}
} // namespace TpcAlignmentLaserRays
