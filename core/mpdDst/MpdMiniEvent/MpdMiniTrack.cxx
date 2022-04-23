//
// MpdMiniTrack holds information about reconstructed track parameters
//

// C++ headers
#include <limits>

// ROOT headers
#include "TMath.h"

// MiniDst headers
#include "MpdMiniMessMgr.h"
#include "MpdMiniTrack.h"

ClassImp(MpdMiniTrack);

//_________________
MpdMiniTrack::MpdMiniTrack()
   : TObject(), fId(0), fChi2(std::numeric_limits<unsigned short>::max()), fPMomentumX(0), fPMomentumY(0),
     fPMomentumZ(0), fGMomentumX(0), fGMomentumY(0), fGMomentumZ(0), fOriginX(0), fOriginY(0), fOriginZ(0),
     fDedx(0),  /* fDedxError(0), */
     fNHits(0), /* fNHitsMax(0), fNHitsDedx(0), */
     fNSigmaPion(std::numeric_limits<short>::min()), fNSigmaKaon(std::numeric_limits<short>::min()),
     fNSigmaProton(std::numeric_limits<short>::min()), fNSigmaElectron(std::numeric_limits<short>::min()), fHitMap(0),
     fBECalClusterIndex(-1), fBTofPidTraitsIndex(-1), fMcTrackId(-1)
{
   // Default constructor
   /* empty */
}

//_________________
MpdMiniTrack::MpdMiniTrack(const MpdMiniTrack &track) : TObject()
{
   // Copy constructor
   fId         = track.fId;
   fChi2       = track.fChi2;
   fPMomentumX = track.fPMomentumX;
   fPMomentumY = track.fPMomentumY;
   fPMomentumZ = track.fPMomentumZ;
   fGMomentumX = track.fGMomentumX;
   fGMomentumY = track.fGMomentumY;
   fGMomentumZ = track.fGMomentumZ;
   fOriginX    = track.fOriginX;
   fOriginY    = track.fOriginY;
   fOriginZ    = track.fOriginZ;
   fDedx       = track.fDedx;
   // fDedxError = track.fDedxError;
   fNHits = track.fNHits;
   // fNHitsMax = track.fNHitsMax;
   // fNHitsDedx = track.fNHitsDedx;
   fNSigmaPion         = track.fNSigmaPion;
   fNSigmaKaon         = track.fNSigmaKaon;
   fNSigmaProton       = track.fNSigmaProton;
   fNSigmaElectron     = track.fNSigmaElectron;
   fHitMap             = track.fHitMap;
   fBECalClusterIndex  = track.fBECalClusterIndex;
   fBTofPidTraitsIndex = track.fBTofPidTraitsIndex;
   fMcTrackId          = track.fMcTrackId;
}

//_________________
MpdMiniTrack::~MpdMiniTrack()
{
   // Destructor
   /* emtpy */
}

//_________________
void MpdMiniTrack::Print(const Char_t *option __attribute__((unused))) const
{
   // Print information about the track
   LOG_INFO << "id: " << id() << " chi2: " << chi2() << "\n"
            << "pMom: " << pMom().X() << " " << pMom().Y() << " " << pMom().Z() << "\n"
            << "gMom: " << gMom().X() << " " << gMom().Y() << " " << gMom().Z() << "\n"
            << "origin: " << origin().X() << " " << origin().Y() << " " << origin().Z() << "\n"
            << "nHits: " << nHits() << "nSigma pi/K/p/e: " << nSigmaPion() << "/" << nSigmaKaon() << "/"
            << nSigmaProton() << "/" << nSigmaElectron() << "\n"
            << "Hit index in ECal/BTof: " << fBECalClusterIndex << "/" << fBTofPidTraitsIndex << "\n"
            << "MC track id: " << fMcTrackId << "\n"
            << endm;
}

//_________________
Float_t MpdMiniTrack::gDCAxy(Float_t x, Float_t y) const
{
   return TMath::Sqrt((fOriginX - x) * (fOriginX - x) + (fOriginY - y) * (fOriginY - y));
}

//_________________
Float_t MpdMiniTrack::gDCA(Float_t x, Float_t y, Float_t z) const
{
   return TMath::Sqrt((fOriginX - x) * (fOriginX - x) + (fOriginY - y) * (fOriginY - y) +
                      (fOriginZ - z) * (fOriginZ - z));
}

//_________________
TVector3 MpdMiniTrack::gDCA(TVector3 pVtx) const
{
   return (origin() - pVtx);
}

//_________________
void MpdMiniTrack::setChi2(Float_t chi2)
{
   fChi2 = ((chi2 * 1000.) > std::numeric_limits<unsigned short>::max() ? std::numeric_limits<unsigned short>::max()
                                                                        : (UShort_t)(TMath::Nint(chi2 * 1000.)));
}

//_________________
void MpdMiniTrack::setDedx(Float_t dEdx)
{
   fDedx = dEdx;
}

//_________________
/*
void MpdMiniTrack::setNHitsMax(Int_t nhits) {
  fNHitsMax = (UChar_t)nhits;
}
*/

//_________________
/*
void MpdMiniTrack::setNHitsPossible(Int_t nhits) {
  // For those who wants to have standard terminology
  setNHitsMax(nhits);
}
*/

//_________________
/*
void MpdMiniTrack::setNHitsDedx(Int_t nhits) {
  fNHitsDedx = (UChar_t)nhits;
}
*/

//_________________
void MpdMiniTrack::setNSigmaPion(Float_t ns)
{
   fNSigmaPion = (fabs(ns * 1000.) > std::numeric_limits<short>::max()
                     ? ((ns > 0) ? std::numeric_limits<short>::max() : std::numeric_limits<short>::min())
                     : (Short_t)(TMath::Nint(ns * 1000.)));
}

//_________________
void MpdMiniTrack::setNSigmaKaon(Float_t ns)
{
   fNSigmaKaon = (fabs(ns * 1000.) > std::numeric_limits<short>::max()
                     ? ((ns > 0) ? std::numeric_limits<short>::max() : std::numeric_limits<short>::min())
                     : (Short_t)(TMath::Nint(ns * 1000.)));
}

//_________________
void MpdMiniTrack::setNSigmaProton(Float_t ns)
{
   fNSigmaProton = (fabs(ns * 1000.) > std::numeric_limits<short>::max()
                       ? ((ns > 0) ? std::numeric_limits<short>::max() : std::numeric_limits<short>::min())
                       : (Short_t)(TMath::Nint(ns * 1000.)));
}

//_________________
void MpdMiniTrack::setNSigmaElectron(Float_t ns)
{
   fNSigmaElectron = (fabs(ns * 1000.) > std::numeric_limits<short>::max()
                         ? ((ns > 0) ? std::numeric_limits<short>::max() : std::numeric_limits<short>::min())
                         : (Short_t)(TMath::Nint(ns * 1000.)));
}

//_________________
TVector3 MpdMiniTrack::gMom(TVector3 pVtx, Float_t const B) const
{
   MpdMiniPhysicalHelix gHelix = helix(B);
   return gHelix.momentumAt(gHelix.pathLength(pVtx), B * kilogauss);
}

//_________________
MpdMiniPhysicalHelix MpdMiniTrack::helix(Float_t const B) const
{
   return MpdMiniPhysicalHelix(gMom(), origin(), B * kilogauss, static_cast<float>(charge()));
}

//_________________
void MpdMiniTrack::setMcTrackIndex(Int_t id)
{
   if (id < 0) {
      fMcTrackId = -1;
   } else {
      fMcTrackId =
         ((id >= 0 && id <= std::numeric_limits<short>::max()) ? (Short_t)id : std::numeric_limits<short>::max());
   }
}
