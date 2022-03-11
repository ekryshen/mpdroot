/**
 * \class MpdMiniPhysicalHelix
 * \brief Helix parametrization for the particle
 *
 * Parametrization of a physical helix that uses ROOT classes
 *
 * \author Grigory Nigmatkulov
 * \email nigmatkulov@gmail.com ; ganigmatkulov@mephi.ru
 * \date May 01, 2020
 */

#ifndef MpdMiniPhysicalHelix_h
#define MpdMiniPhysicalHelix_h

// ROOT headers
#include "TVector3.h"

// MiniDst headers
#include "MpdMiniHelix.h"

//_________________
class MpdMiniPhysicalHelix : public MpdMiniHelix {

public:
   /// Empty constructor
   MpdMiniPhysicalHelix();
   /// Constructor with momentum, origin, signed Magnetic Field
   /// and Charge of particle (+/- 1)
   MpdMiniPhysicalHelix(const TVector3 &, const TVector3 &, Double_t, Double_t);
   /// Constructor with Curvature, dip angle, phase, origin, h
   MpdMiniPhysicalHelix(Double_t, Double_t, Double_t, const TVector3 &, Int_t h = -1);
   /// Destructor
   ~MpdMiniPhysicalHelix();

   /// Return the momentum at origin
   ///    \param bField magnetic field
   TVector3 momentum(Double_t) const;
   /// Return momemtum at S
   TVector3 momentumAt(Double_t, Double_t) const;
   /// Return charge of a particle
   Int_t charge(Double_t) const;
   /// 2d DCA to x,y point signed relative to curvature
   Double_t curvatureSignedDistance(Double_t x, Double_t y);
   /// 2d DCA to x,y point signed relative to rotation
   Double_t geometricSignedDistance(Double_t x, Double_t y);
   /// 3d DCA to 3d point signed relative to curvature
   Double_t curvatureSignedDistance(const TVector3 &);
   /// 3d DCA to 3d point signed relative to rotation
   Double_t geometricSignedDistance(const TVector3 &);

   ClassDef(MpdMiniPhysicalHelix, 1)
};

#endif // #define MpdMiniPhysicalHelix_h
