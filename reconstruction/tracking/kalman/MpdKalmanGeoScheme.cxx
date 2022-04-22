/// \class MpdKalmanGeoScheme
///
/// Geometry manager for MPD Kalman filter
/// \author Alexander Zinchenko (LHEP, JINR, Dubna)

#include "MpdKalmanGeoScheme.h"
#include "MpdKalmanHit.h"
//#include "FairRootManager.h"
//#include <TMath.h>
#include <Riostream.h>
#include <map>

using namespace std;

//__________________________________________________________________________
MpdKalmanGeoScheme::MpdKalmanGeoScheme() : TObject()
{
   /// Default constructor
}

//__________________________________________________________________________
MpdKalmanGeoScheme::~MpdKalmanGeoScheme()
{
   /// Destructor
   // FairRootManager *manager = FairRootManager::Instance();
   // manager->Write();
   fPosMap.clear();
   fNormMap.clear();
   fSizeMap.clear();
   fDetMap.clear();
}

//__________________________________________________________________________
Bool_t MpdKalmanGeoScheme::Exist(Int_t detID)
{
   // Check if detector detID exists in geo scheme

   if (fPosMap.find(detID) == fPosMap.end()) return kFALSE;
   return kTRUE;
}

//__________________________________________________________________________
TVector3 MpdKalmanGeoScheme::GlobalPos(const MpdKalmanHit *hit)
{
   // Get global coordinates of the hit

   return fPosMap[hit->GetDetectorID()];
}

//__________________________________________________________________________
TVector3 MpdKalmanGeoScheme::Normal(const MpdKalmanHit *hit)
{
   // Get normal direction to the plane hit

   // return TVector3(0.,0.,0.); // just dummy
   return fNormMap[hit->GetDetectorID()];
}

//__________________________________________________________________________
TVector2 MpdKalmanGeoScheme::Size(const MpdKalmanHit *hit)
{
   // Get normal direction to the plane hit

   return fSizeMap[hit->GetDetectorID()];
}

//__________________________________________________________________________
Int_t MpdKalmanGeoScheme::DetId(TString detName)
{
   // Get detector ID from its name

   if (fDetMap.find(detName) == fDetMap.end()) return -1;
   return fDetMap[detName];
}

//__________________________________________________________________________
TString MpdKalmanGeoScheme::Path(Int_t detId)
{
   // Get detector path from its ID

   if (fPathMap.find(detId) == fPathMap.end()) return -1;
   return fPathMap[detId];
}

//__________________________________________________________________________
void MpdKalmanGeoScheme::SetGlobalPos(const MpdKalmanHit *hit, TVector3 pos, Bool_t erase)
{
   // Put hit in GeoScheme

   Int_t detID = hit->GetDetectorID();
   SetGlobalPos(detID, pos, erase);
}

//__________________________________________________________________________
void MpdKalmanGeoScheme::SetGlobalPos(Int_t detID, TVector3 pos, Bool_t erase)
{
   // Put detector position in GeoScheme

   if (erase) fPosMap.erase(detID);
   fPosMap.insert(pair<Int_t, TVector3>(detID, pos));
}

//__________________________________________________________________________
void MpdKalmanGeoScheme::SetNormal(Int_t detID, TVector3 norm, Bool_t erase)
{
   // Put detector normal in GeoScheme

   if (erase) fNormMap.erase(detID);
   fNormMap.insert(pair<Int_t, TVector3>(detID, norm));
}

//__________________________________________________________________________
void MpdKalmanGeoScheme::SetSize(Int_t detID, TVector2 size, Bool_t erase)
{
   // Put detector dimensions in GeoScheme

   if (erase) fSizeMap.erase(detID);
   fSizeMap.insert(pair<Int_t, TVector2>(detID, size));
}

//__________________________________________________________________________
void MpdKalmanGeoScheme::SetDetId(TString detName, Int_t detID)
{
   // Store detector ID with name detName

   fDetMap.insert(pair<TString, Int_t>(detName, detID));
}

//__________________________________________________________________________
void MpdKalmanGeoScheme::SetPath(Int_t detID, TString path, Bool_t erase)
{
   // Store detector path with detID

   if (erase) fPathMap.erase(detID);
   fPathMap.insert(pair<Int_t, TString>(detID, path));
}

//__________________________________________________________________________
void MpdKalmanGeoScheme::Reset()
{
   // Reset maps

   fPosMap.clear();
   fNormMap.clear();
   fSizeMap.clear();
   fDetMap.clear();
   fPathMap.clear();
}

ClassImp(MpdKalmanGeoScheme);
