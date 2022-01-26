//--------------------------------------------------------------------
//
// Description:
//      MPD EMC positions of track extrapolation to ECAL
//
//
// Author List:
//      D.Peresunko
//
//--------------------------------------------------------------------

#ifndef MPDEMCTRACKEXTRAP_H
#define MPDEMCTRACKEXTRAP_H
#include "TObject.h"
class MpdEmcClusterKI;
class MpdEmcTrackExtrap: public TObject
{
 public:
  /** Standard constructor**/
  MpdEmcTrackExtrap(int trckInd, double xTop, double yTop, double zTop, double xD, double yD, double zD);

  /** Destructor **/
  ~MpdEmcTrackExtrap() {}

  double Distance(MpdEmcClusterKI* clu) const;

  void DistanceDphiDz(MpdEmcClusterKI* clu, double& dphi, double& dz) const;

  void DistanceLongPerp(MpdEmcClusterKI* clu, double& dlong, double& dPerp) const;

  int GetTrackIndex() const { return fIndex; }

 protected:
  int fIndex; // Track index
  double fX1; // Extrapolated to ECAL surface x
  double fY1; // Extrapolated to ECAL surface y
  double fZ1; // Extrapolated to ECAL surface z
  double fX2; // Extrapolated to 6X0 depth x
  double fY2; // Extrapolated to 6X0 depth y
  double fZ2; // Extrapolated to 6X0 depth z

  ClassDef(MpdEmcTrackExtrap, 1);
};

#endif
