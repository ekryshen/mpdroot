#ifndef MPDKALMANGEOSCHEME_H
#define MPDKALMANGEOSCHEME_H

/// \ingroup rec
/// \class MpdKalmanGeoScheme
/// \brief Geometry manager for MPD Kalman filter
///
/// \author Alexander Zinchenko, LHEP JINR Dubna

class MpdKalmanHit;
#include <TString.h>
#include <TVector2.h>
#include <TVector3.h>
#include <map>

class MpdKalmanGeoScheme : public TObject {

public:
   MpdKalmanGeoScheme();          ///< Default ctor
   virtual ~MpdKalmanGeoScheme(); ///< Destructor

   Bool_t   Exist(Int_t detID);                 ///< check if detector detID exists in geo scheme
   TVector3 GlobalPos(const MpdKalmanHit *hit); ///< get global coordinates of the hit
   TVector3 Normal(const MpdKalmanHit *hit);    ///< get normal direction to the plane hit
   TVector2 Size(const MpdKalmanHit *hit);      ///< get detector dimensions
   Int_t    DetId(TString detName);             ///< get detector ID from its name
   TString  Path(Int_t detId);                  ///< get detector path from its ID
   void     SetGlobalPos(const MpdKalmanHit *hit, TVector3 pos, Bool_t erase = kFALSE); ///< put hit in GeoScheme
   void     SetGlobalPos(Int_t detID, TVector3 pos, Bool_t erase = kFALSE); ///< put detector position in GeoScheme
   void     SetNormal(Int_t detID, TVector3 norm, Bool_t erase = kFALSE);   ///< put detector normal in GeoScheme
   void     SetSize(Int_t detID, TVector2 size, Bool_t erase = kFALSE);     ///< put detector dimensions in GeoScheme
   void     SetDetId(TString detName, Int_t detID);                         ///< store detector ID with name detName
   void     SetPath(Int_t detID, TString path, Bool_t erase = kFALSE);      ///< store detector path with detID
   void     Reset();                                                        ///< reset maps

protected:
   // virtual void Finish();

private:
   std::map<Int_t, TVector3> fPosMap;
   std::map<Int_t, TVector3> fNormMap;
   std::map<Int_t, TVector2> fSizeMap;
   std::map<TString, Int_t>  fDetMap;
   std::map<Int_t, TString>  fPathMap;

   ClassDef(MpdKalmanGeoScheme, 1);
};
#endif
