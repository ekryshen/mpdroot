//
// MpdMiniBECalCluster stores cluster information in ECal
//

// C++ headers
#include <limits>

// ROOT headers
#include "TString.h"

// MiniDst headers
#include "MpdMiniMessMgr.h"
#include "MpdMiniBECalCluster.h"

ClassImp(MpdMiniBECalCluster);

//_________________
MpdMiniBECalCluster::MpdMiniBECalCluster()
   : TObject(), fEnergy(0), fEcore(0), fEcore1p(0), fEcore2p(0), fTime(0), fX(0), fY(0), fZ(0), fdPhi(0), fdZ(0),
     fTrackId(-1), /* fDisp(0), */ fLambda1(0), fLambda2(0), fChi2(0), fNExLM(0)
{
   // Default constructor
   if (!fDigitIds.empty()) fDigitIds.clear();
   if (!fDigitEDeps.empty()) fDigitEDeps.clear();
   if (!fMcTrackIds.empty()) fMcTrackIds.clear();
   if (!fMcTrackEDeps.empty()) fMcTrackEDeps.clear();
}

//_________________
MpdMiniBECalCluster::MpdMiniBECalCluster(const MpdMiniBECalCluster &c)
   : TObject(), fDigitIds(c.fDigitIds), fDigitEDeps(c.fDigitEDeps), fMcTrackIds(c.fMcTrackIds),
     fMcTrackEDeps(c.fMcTrackEDeps), fEnergy(c.fEnergy), fEcore(c.fEcore), fEcore1p(c.fEcore1p), fEcore2p(c.fEcore2p),
     fTime(c.fTime), fX(c.fX), fY(c.fY), fZ(c.fZ), fdPhi(c.fdPhi), fdZ(c.fdZ), fTrackId(c.fTrackId),
     /* fDisp(c.fDisp), */ fLambda1(c.fLambda1), fLambda2(c.fLambda2), fChi2(c.fChi2), fNExLM(c.fNExLM)
{
   // Copy constructor
   /* empty */
}

//_________________
MpdMiniBECalCluster::~MpdMiniBECalCluster()
{
   // Destructor
   if (!fDigitIds.empty()) fDigitIds.clear();
   if (!fDigitEDeps.empty()) fDigitEDeps.clear();
   if (!fMcTrackIds.empty()) fMcTrackIds.clear();
   if (!fMcTrackEDeps.empty()) fMcTrackEDeps.clear();
}

//__________________
void MpdMiniBECalCluster::Print(const Char_t *option __attribute__((unused))) const
{
   LOG_INFO << "nDigits: " << fDigitIds.size() << " nMcTracks: " << fMcTrackIds.size() << " Energy: " << fEnergy
            << " Ecore: " << fEcore << " Ecore1p: " << fEcore1p << " Ecore2p: " << fEcore2p << " Time: " << fTime
            << " x/y/z: " << Form("%4.2f/%4.2f/%4.2f", fX, fY, fZ) << " dPhi/dZ: " << Form("%4.2f/%4.2f", fdPhi, fdZ)
            << " TrackId: "
            << fTrackId
            /* << " Disp: " << fDisp */
            << " Lam1/Lam2: " << Form("%4.2f/%4.2f", fLambda1, fLambda2) << " chi2: " << chi2() << " NLM: " << nlm()
            << endm;
}

//_________________
void MpdMiniBECalCluster::digitParams(UShort_t i, Int_t &digId, Float_t &eneDep)
{
   digId  = (i < fDigitIds.size()) ? fDigitIds.at(i) : -1;
   eneDep = (i < fDigitEDeps.size()) ? fDigitEDeps.at(i) : -999.;
}

//_________________
void MpdMiniBECalCluster::mcTrackParams(UShort_t i, Int_t &mcTrkId, Float_t &eneDep)
{
   mcTrkId = (i < fMcTrackIds.size()) ? fMcTrackIds.at(i) : -1;
   eneDep  = (i < fMcTrackEDeps.size()) ? fMcTrackEDeps.at(i) : -999.;
}

//_________________
void MpdMiniBECalCluster::setDigitIds(std::vector<Int_t> ids)
{
   for (UInt_t i = 0; i < ids.size(); i++) {
      fDigitIds.push_back((ids.at(i) < 0 || ids.at(i) > std::numeric_limits<UShort_t>::max())
                             ? (UShort_t)ids.at(i)
                             : std::numeric_limits<UShort_t>::max());
   } // for (UInt_t i=0; i<ids.size(); i++)
}

//_________________
void MpdMiniBECalCluster::setDigitEDeps(std::vector<Float_t> eDeps)
{
   for (UInt_t i = 0; i < eDeps.size(); i++) {
      fDigitEDeps.push_back(eDeps.at(i));
   } // for (UInt_t i=0; i<ids.size(); i++)
}

//_________________
void MpdMiniBECalCluster::addDigit(Int_t id, Float_t eDep)
{
   fDigitIds.push_back((id < 0 || id > std::numeric_limits<UShort_t>::max()) ? (UShort_t)id
                                                                             : std::numeric_limits<UShort_t>::max());
   fDigitEDeps.push_back(eDep);
}

//_________________
void MpdMiniBECalCluster::setMcTrackIds(std::vector<Int_t> ids)
{
   for (UInt_t i = 0; i < ids.size(); i++) {
      fMcTrackIds.push_back((ids.at(i) < 0 || ids.at(i) > std::numeric_limits<UShort_t>::max())
                               ? (UShort_t)ids.at(i)
                               : std::numeric_limits<UShort_t>::max());
   } // for (UInt_t i=0; i<ids.size(); i++)
}

//_________________
void MpdMiniBECalCluster::setMcTrackEDeps(std::vector<Float_t> eDeps)
{
   for (UInt_t i = 0; i < eDeps.size(); i++) {
      fMcTrackEDeps.push_back(eDeps.at(i));
   } // for (UInt_t i=0; i<eDeps.size(); i++)
}

//_________________
void MpdMiniBECalCluster::addMcTrack(Int_t id, Float_t eDep)
{
   fMcTrackIds.push_back((id < 0 || id > std::numeric_limits<UShort_t>::max()) ? (UShort_t)id
                                                                               : std::numeric_limits<UShort_t>::max());
   fMcTrackEDeps.push_back(eDep);
}

//_________________
void MpdMiniBECalCluster::setTrackId(Int_t id)
{
   if (id < std::numeric_limits<Short_t>::min() || id > std::numeric_limits<Short_t>::max()) {
      fTrackId = std::numeric_limits<Short_t>::min();
   } else {
      fTrackId = (Short_t)id;
   }
}

//_________________
void MpdMiniBECalCluster::setChi2(Float_t chi2)
{
   if ((chi2 * 100.f) < 0 || (chi2 * 100.) > std::numeric_limits<UShort_t>::max()) {
      fChi2 = std::numeric_limits<UShort_t>::max();
   } else {
      fChi2 = (UShort_t)(chi2 * 100);
   }
}

//_________________
void MpdMiniBECalCluster::setNLM(Int_t nlm)
{
   if (nlm < 0) {
      fNExLM = 0;
   } else if (nlm > std::numeric_limits<UChar_t>::max()) {
      fNExLM = std::numeric_limits<UChar_t>::max();
   } else {
      fNExLM = (UChar_t)nlm;
   }
}
