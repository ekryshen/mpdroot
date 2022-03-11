/**
 * \class MpdMiniTrack
 * \brief Holds information about reconstructed track parameters
 *
 * The class stores information about tracks reconstructed in TPC
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \email nigmatkulov@gmail.com ; ganigmatkulov@mephi.ru
 * \date May 01, 2020
 */

#ifndef MpdMiniTrack_h
#define MpdMiniTrack_h

// C++ headers
#include <cmath>

// ROOT headers
#include "TObject.h"
#include "TVector3.h"
#include "TBits.h"

// MiniDst headers
#include "MpdMiniHelix.h"
#include "MpdMiniPhysicalHelix.h"
#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

//_________________
class MpdMiniTrack : public TObject {

public:
   /// Default constructor
   MpdMiniTrack();
   /// Copy constructor
   MpdMiniTrack(const MpdMiniTrack &track);
   /// Destructor
   virtual ~MpdMiniTrack();
   /// Print track parameters
   virtual void Print(const Char_t *option = "") const;

   //
   // Getters
   //

   /// Return unique Id of the track
   Int_t id() const { return fId; }
   /// Return chi2 of the track
   Float_t chi2() const { return fChi2 / 1000.f; }
   /// Return momentum (GeV/c) of the primary track. Return (0.,0.,0.) if not primary track
   TVector3 pMom() const { return TVector3(fPMomentumX, fPMomentumY, fPMomentumZ); }
   /// Return momentum (GeV/c) of the global tracks at the point of DCA to the primary vertex
   TVector3 gMom() const { return TVector3(fGMomentumX, fGMomentumY, fGMomentumZ); }
   /// Return space coordinate (cm) of DCA to the primary vertex
   TVector3 origin() const { return TVector3(fOriginX, fOriginY, fOriginZ); }
   /// Return tranverse momentum (GeV/c) of the global track
   Float_t gPt() const { return gMom().Perp(); }
   /// Return total momentum (GeV/c) of the global track
   Float_t gPtot() const { return gMom().Mag(); }
   /// Return transverse momentum (GeV/c) of the primary track
   Float_t pPt() const { return isPrimary() ? pMom().Perp() : 0.; }
   /// Return total momentum (GeV/c) of the primary track
   Float_t pPtot() const { return isPrimary() ? pMom().Mag() : 0.; }
   /// Global momentum at point of DCA to pVtx, B should be in kilogauss
   TVector3 gMom(TVector3 pVtx, Float_t B) const;
   /// Helix at point of DCA to MpdMiniEvent::mPrimaryVertex
   MpdMiniPhysicalHelix helix(Float_t const B) const;

   /// Return signed distance in x direction (cm) between the value
   /// and x position of the DCA point (gDCAx - x)
   Float_t gDCAx(Float_t pVtxX) const { return (fOriginX - pVtxX); }
   /// Return signed distance in y direction (cm) between the value
   /// and y position of the DCA point (gDCAy - y)
   Float_t gDCAy(Float_t pVtxY) const { return (fOriginY - pVtxY); }
   /// Return signed distance in z direction (cm) between the value
   /// and z position of the DCA point (gDCAz - z)
   Float_t gDCAz(Float_t pVtxZ) const { return (fOriginZ - pVtxZ); }
   /// Return distance in xy direction (cm) between the (x,y) point
   /// and the DCA point to primary vertex
   Float_t gDCAxy(Float_t pVtxX, Float_t pVtxY) const;
   /// Return distance in xyz direction (cm) between the (x,y,z) point
   /// and the DCA point to primary vertex
   Float_t gDCA(Float_t pVtxX, Float_t pVtxY, Float_t pVtxZ) const;
   /// Return 3-vector (distance) between the DCA point and the point (gDCA - point)
   TVector3 gDCA(TVector3 pVtx) const;
   /// Return charge of the track (encoded in nHitsFit as: nHitsFit * charge)
   Int_t charge() const { return (fNHits > 0) ? 1 : -1; }
   /// Return number of hits
   Int_t nHits() const { return (fNHits > 0) ? (Int_t)fNHits : (Int_t)(-1 * fNHits); }

   // Return number of hits possible
   // Int_t   nHitsMax() const               { return (Int_t)fNHitsMax; }

   // Return number of hits possible
   // Int_t   nHitsPoss() const              { return nHitsMax(); }

   // Return number of hits used for dE/dx measurement
   // Int_t   nHitsDedx() const              { return (Int_t)fNHitsDedx; }

   /// Return dE/dx (GeV/cm) of the track
   Float_t dEdx() const { return fDedx; }

   // Return dE/dx error of the track
   // Float_t dEdxError() const              { return fDedxError; }

   /// Return nSigma(pion)
   Float_t nSigmaPion() const { return (Float_t)fNSigmaPion / 1000.f; }
   /// Return nSigma(kaon)
   Float_t nSigmaKaon() const { return (Float_t)fNSigmaKaon / 1000.f; }
   /// Return nSigma(proton)
   Float_t nSigmaProton() const { return (Float_t)fNSigmaProton / 1000.f; }
   /// Return nSigma(electron)
   Float_t nSigmaElectron() const { return (Float_t)fNSigmaElectron / 1000.f; }

   /// Return hit map
   ULong64_t hitMap() const { return fHitMap; }

   /// Return if track has TOF hit
   Bool_t isBTofTrack() const { return (fBTofPidTraitsIndex < 0) ? kFALSE : kTRUE; }
   /// Return if track was reconstructed in ECal cluster
   Bool_t isBECalTrack() const { return (fBECalClusterIndex < 0) ? kFALSE : kTRUE; }

   /// Return if track is primary
   Bool_t isPrimary() const { return (pMom().Mag() > 0); }

   /// Return index of the associated barrel ECal cluster (-1 if not match)
   Int_t bECalClusterIndex() const { return fBECalClusterIndex; }
   /// Return index of the associated barrel TofPidTrait (-1 if not match)
   Int_t bTofPidTraitsIndex() const { return fBTofPidTraitsIndex; }

   /// Return index of MpdMiniMcTrack corresponding
   /// to current MpdMiniTrack (-1 if not match)
   Int_t mcTrackIndex() const { return (Int_t)fMcTrackId; }
   /// Check if miniTrack has corresponding miniMcTrack
   Bool_t hasMcTrack() const { return (fMcTrackId < 0) ? kFALSE : kTRUE; }

   //
   // Setters
   //

   /// Set track ID
   void setId(Int_t id) { fId = (UShort_t)id; }
   /// Set chi2 of the track
   void setChi2(Float_t chi2);
   /// Set momentum of the primary track
   void setPrimaryMomentum(Double_t px, Double_t py, Double_t pz)
   {
      fPMomentumX = (Float_t)px;
      fPMomentumY = (Float_t)py;
      fPMomentumZ = (Float_t)pz;
   }
   /// Set momentum of the primary track
   void setPrimaryMomentum(Float_t px, Float_t py, Float_t pz)
   {
      fPMomentumX = px;
      fPMomentumY = py;
      fPMomentumZ = pz;
   }
   /// Set momentum of the primary track
   void setPrimaryMomentum(TVector3 mom)
   {
      fPMomentumX = (Float_t)mom.X();
      fPMomentumY = (Float_t)mom.Y();
      fPMomentumZ = (Float_t)mom.Z();
   }
   /// Set momentum of the global track
   void setGlobalMomentum(Double_t px, Double_t py, Double_t pz)
   {
      fGMomentumX = (Float_t)px;
      fGMomentumY = (Float_t)py;
      fGMomentumZ = (Float_t)pz;
   }
   /// Set momentum of the global track
   void setGlobalMomentum(Float_t px, Float_t py, Float_t pz)
   {
      fGMomentumX = px;
      fGMomentumY = py;
      fGMomentumZ = pz;
   }
   /// Set momentum of the global track
   void setGlobalMomentum(TVector3 mom)
   {
      fGMomentumX = (Float_t)mom.X();
      fGMomentumY = (Float_t)mom.Y();
      fGMomentumZ = (Float_t)mom.Z();
   }
   /// Set origin of the track (DCA point to the primary vertex)
   void setOrigin(Double_t x, Double_t y, Double_t z)
   {
      fOriginX = (Float_t)x;
      fOriginY = (Float_t)y;
      fOriginZ = (Float_t)z;
   }
   /// Set origin of the track (DCA point to the primary vertex)
   void setOrigin(Float_t x, Float_t y, Float_t z)
   {
      fOriginX = x;
      fOriginY = y;
      fOriginZ = z;
   }
   /// Set origin of the track (DCA point to the primary vertex)
   void setOrigin(TVector3 orig)
   {
      fOriginX = (Float_t)orig.X();
      fOriginY = (Float_t)orig.Y();
      fOriginZ = (Float_t)orig.Z();
   }

   /// Set dE/dx of the track
   void setDedx(Float_t dEdx);

   // Set dE/dx error of the track
   // void setDedxError(Float_t dEdxError)     { fDedxError = dEdxError; }

   /// Set nHitsFit ( charge * nHitsFit )
   void setNHits(Int_t nhits) { fNHits = (Char_t)nhits; }

   // Set nHitsPoss
   // void setNHitsPossible(Int_t nhits);

   // Set nHitsPoss
   // void setNHitsMax(Int_t nhits);

   /// Set nHitsDedx
   void setNHitsDedx(Int_t nhits);
   /// Set nSigma(pion)
   void setNSigmaPion(Float_t ns);
   /// Set nSigma(kaon)
   void setNSigmaKaon(Float_t ns);
   /// Set nSigma(proton)
   void setNSigmaProton(Float_t ns);
   /// Set nSigma(electron)
   void setNSigmaElectron(Float_t ns);
   /// Set hit map
   void setHitMap(ULong64_t map) { fHitMap = map; }

   /// Set index of the corresponding miniMcTrack
   void setMcTrackIndex(Int_t index);

   /// Set index to ECal cluster
   void setBECalClusterIndex(Int_t index) { fBECalClusterIndex = (Short_t)index; }
   /// Set index to BTOF PID traits
   void setBTofPidTraitsIndex(Int_t index) { fBTofPidTraitsIndex = (Short_t)index; }

protected:
   /// Unique track ID
   UShort_t fId;
   /// Chi2 of the track (encoding = chi2*1000)
   UShort_t fChi2;
   /// Px momentum (GeV/c) of the primary track ( 0 if not primary )
   Float_t fPMomentumX;
   /// Py momentum (GeV/c) of the primary track ( 0 if not primary )
   Float_t fPMomentumY;
   /// Pz momentum (GeV/c) of the primary track ( 0 if not primary )
   Float_t fPMomentumZ;
   /// Px component of the momentum (GeV/c) of the global track
   /// at DCA to primary vertex
   Float_t fGMomentumX;
   /// Py component of the momentum (GeV/c) of the global track
   /// at DCA to primary vertex
   Float_t fGMomentumY;
   /// Pz component of the momentum (GeV/c) of the global track
   /// at DCA to primary vertex
   Float_t fGMomentumZ;
   /// Track origin x in cm (at DCAx to the primary vertex)
   Float_t fOriginX;
   /// Track origin y in cm (at DCAy to the primary vertex)
   Float_t fOriginY;
   /// Track origin z in cm (DCAy to the primary vertex)
   Float_t fOriginZ;

   /// dE/dx in arbitrary units
   Float16_t fDedx;

   // dE/dx error in arbitrary units
   // Float16_t fDedxError;

   /// Charge * nHits
   Char_t fNHits;

   // Possible number of hits (in TPC)
   // UChar_t  fNHitsMax;

   // Number of hits used for dE/dx estimation (in TPC)
   // UChar_t  fNHitsDedx;

   /// nSigma(pion)  (encoding = nsigma * 1000)
   Short_t fNSigmaPion;
   /// nSigma(kaon)  (encoding = nsigma * 1000)
   Short_t fNSigmaKaon;
   /// nSigma(proton)  (encoding = nsigma * 1000)
   Short_t fNSigmaProton;
   /// nSigma(electron)  (encoding = nsigma * 1000)
   Short_t fNSigmaElectron;
   /// Track hit map
   ULong64_t fHitMap;

   /// Index of the barrel ECal cluster in the event (-1 if no match)
   Short_t fBECalClusterIndex;
   /// Index of the BTOF pidTratis in the event (-1 if no match)
   Short_t fBTofPidTraitsIndex;
   /// Index of miniMcTrack that correponds to miniTrack (aka IdTruth).
   /// (-1 if no match)
   Short_t fMcTrackId;

   ClassDef(MpdMiniTrack, 3)
};

#endif // #define MpdMiniTrack_h
