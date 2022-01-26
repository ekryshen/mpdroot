#include "MpdEmcTrackExtrap.h"
#include "MpdEmcClusterKI.h"

MpdEmcTrackExtrap::MpdEmcTrackExtrap(int trckInd, double xTop, double yTop, double zTop, double xD, double yD,
                                     double zD)
  : fIndex(trckInd), fX1(xTop), fY1(yTop), fZ1(zTop), fX2(xD), fY2(yD), fZ2(zD)
{
}

double MpdEmcTrackExtrap::Distance(MpdEmcClusterKI* clu) const
{
  double dx = fX2 - clu->GetX(), dy = fY2 - clu->GetY(), dz = fZ2 - clu->GetZ();
  return sqrt(dx * dx + dy * dy + dz * dz);
}

void MpdEmcTrackExtrap::DistanceDphiDz(MpdEmcClusterKI* clu, double& dphi, double& dz) const
{
  double dx = fX2 - clu->GetX(), dy = fY2 - clu->GetY();
  dphi = (dx * fY2 - dy * fX2) / sqrt(fX2 * fX2 + fY2 * fY2); // TODO: Check sign!
  dz = fZ2 - clu->GetZ();
}

void MpdEmcTrackExtrap::DistanceLongPerp(MpdEmcClusterKI* clu, double& dLong, double& dPerp) const
{
  double dx = fX2 - clu->GetX(), dy = fY2 - clu->GetY(), dz = fZ2 - clu->GetZ();
  double dxT = fX2 - fX1, dyT = fY2 - fY1, dzT = fZ2 - fZ1;
  double dl = sqrt(dxT * dxT + dyT * dyT + dzT * dzT);
  if (dl > 0) {
    dLong = (dx * dxT + dy * dyT + dz * dzT) / dl; // scalar product
    double xi = dy * dzT - dz * dyT;
    double yi = dz * dxT - dz * dzT;
    dPerp = TMath::Sign(sqrt(xi * xi + yi * yi) / dl, xi * fX2 + yi * fY2); // vector product
  }
}

ClassImp(MpdEmcTrackExtrap)

