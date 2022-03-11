/**
 * Hold information about BTOF-matched tracks
 *
 * The class stores information about tracks that matched
 * the Barrel Time-of-Flight detector
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI), Pavel Batyuk (JINR)
 * \email ganigmatkulov@mephi.ru ; nigmatkulov@gmail.com ; pavel.batyuk@jinr.ru
 * \date May 01, 2020
 */

#ifndef MpdMiniBTofPidTraits_h
#define MpdMiniBTofPidTraits_h

// C++ headers
#include <limits>

// ROOT headers
#include "TObject.h"
#include "TVector3.h"

//_________________
class MpdMiniBTofPidTraits : public TObject {

public:
   /// Default constructor
   MpdMiniBTofPidTraits();
   /// Copy constructor
   MpdMiniBTofPidTraits(const MpdMiniBTofPidTraits &traits);
   /// Destructor
   virtual ~MpdMiniBTofPidTraits();
   /// Print TOF PID traits information
   virtual void Print(const Char_t *option = "") const;

   //
   // Getters
   //

   /// Return index of the associated track
   Int_t trackIndex() const { return fTrackIndex; }
   /// Return index of the associated hit
   Int_t hitIndex() const { return fHitIndex; }
   /// Return beta (compression = beta * 20000)
   Float_t beta() const { return (Float_t)fBTofBeta / 20000.f; }
   /// Return track momentum
   TVector3 p() const { return TVector3(fPx, fPy, fPz); }
   /// Return mass square
   Float_t massSqr() const;

   //
   // Setters
   //

   /// Set associated track index
   void setTrackIndex(Int_t idx) { fTrackIndex = (idx > std::numeric_limits<Short_t>::max()) ? -1 : (Short_t)idx; }
   /// Set associated tof hit index
   void setHitIndex(Int_t idx) { fHitIndex = (idx > std::numeric_limits<Short_t>::max()) ? -1 : (Short_t)idx; }
   /// Set beta
   void setBeta(Float_t beta);
   /// Set momentum
   void setMomentum(TVector3 p)
   {
      fPx = p.X();
      fPy = p.Y();
      fPz = p.Z();
   }

private:
   /// Index to the associated track in the event (-1 if no matching)
   Short_t fTrackIndex;
   /// Index to the associated hit in the event (-1 if no matching)
   Short_t fHitIndex;
   /// Beta * 20000
   UShort_t fBTofBeta;
   /// Px of the track (GeV/c)
   Float16_t fPx;
   /// Py of the track (GeV/c)
   Float16_t fPy;
   /// Pz of the track (GeV/c)
   Float16_t fPz;

   ClassDef(MpdMiniBTofPidTraits, 3);
};

#endif // #define MpdMiniBTofPidTraits_h
