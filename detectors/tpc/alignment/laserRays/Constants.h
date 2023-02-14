#ifndef CONSTANTS_HH
#define CONSTANTS_HH
#include <cmath>

namespace TpcAlignmentLaserRays {
const double PI(::std::atan(1.) * 4.);
const double rPI(6. / PI);
#define M_ 4
#define N_ 3
#define L_ 3 * M_ // 4*3 = 12
#define LL L_ *L_ // 12*12 = 144
#define SECTORS 24
// vCorrectionMatrix[24*L_] = vCorrectionMatrix[288]
// M[24*LL] = M[3456]
// R[24*L_] = R[288]

} // namespace TpcAlignmentLaserRays

#endif