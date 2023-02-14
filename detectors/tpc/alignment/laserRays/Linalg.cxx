#include "Linalg.h"
#include <cmath>
#include <stdexcept>
#include <string>
#include <climits>

namespace TpcAlignmentLaserRays {

void Linalg::lsolve(std::vector<double> &vA, std::vector<double> &B, std::vector<double> &X)
{
   double r{0.};
   size_t m{0};
   size_t n{B.size()};
   double delta{1.E-10};
   for (size_t i = 0; i < n; i++) {
      m = 0;
      r = 0.;
      for (size_t j = 0; j < n; j++) {
         if (fabs(vA[i + j * n]) > r) {
            m = j;
            r = fabs(vA[i + j * n]);
         }
      }
      r = vA[i + m * n];
      if (fabs(r) < delta) { // r~=0.
         std::string err("R is close to zero: ");
         err.append(std::to_string(r)).append(" < ").append(std::to_string(delta));
         throw std::runtime_error(err.c_str());
      }
      for (size_t j = 0; j < n; j++) {
         vA[i + j * n] = vA[i + j * n] / r;
      }
      B[i]          = B[i] / r;
      vA[i + m * n] = 1.;
      for (size_t k = i + 1; k < n; k++) {
         r = vA[k + m * n];
         for (int j = 0; j < n; j++) {
            vA[k + j * n] = vA[k + j * n] - r * vA[i + j * n];
         }
         B[k]          = B[k] - r * B[i];
         vA[k + m * n] = 0.;
      }
   } // end of first for-loop

   if (n > 0) {
      if (n > INT_MAX) {
         throw std::overflow_error(std::string("size of B-vector is quite big: ").append(std::to_string(n)).c_str());
      }
      int intN = static_cast<int>(n);
      for (int i = intN - 1; i >= 0; --i) {
         m = 0;
         for (size_t j = 0; j < n; j++) {
            if (fabs(vA[i + j * n] - 1.) < delta) { // vA[i + j * n] ~= 1.
               m = j;
            }
         }
         X[m] = B[i];
         if (i > 0) {
            for (int k = i - 1; i > 0 && k >= 0; --k) {
               r             = vA[k + m * n];
               B[k]          = B[k] - r * B[i];
               vA[k + m * n] = 0.;
            }
         }
      }
   }
}
} // namespace TpcAlignmentLaserRays
