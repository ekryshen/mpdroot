//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Found Clusters
//
//
//
// Environment:
//      Software developed for MPD at NICA
//
// Author List:
//      Salmin Roman          (original author)
//
//-----------------------------------------------------------

#ifndef TPCFOUNDCLUSTER_HH
#define TPCFOUNDCLUSTER_HH

// Base Class Headers ----------------
#include "TObject.h"
#include <TVector2.h>
#include <TVector3.h>

#include <TpcOriginated.h>

class TpcFoundCluster : public TObject, public TpcOriginated {
private:
   TVector3 _localPos; // position in Sector reference frame
   Int_t    _sectorID; // Id of sector
   UInt_t   _MCTrackID;

   Double_t _charge; // cluster charge, added by Sergey Merts

   UInt_t   _peak;      // value of peak for checking
   Double_t _meanValey; // mean value of valey for checking
public:
   TpcFoundCluster();

   TpcFoundCluster(const TVector3 &localPos, Int_t sectorID, UInt_t peak, Double_t valey, Double_t charge,
                   UInt_t MCTrackID); // charge added by Sergey Merts
   virtual ~TpcFoundCluster();

   const TVector3 &localPos() const { return _localPos; }
   Int_t           sectorID() const { return _sectorID; }
   Double_t        meanValey() const { return _meanValey; }
   UInt_t          peak() const { return _peak; }
   Double_t        charge() const { return _charge; }
   UInt_t          MCTrackID() const { return _MCTrackID; }

   TVector3 globalPos() const;

   ClassDef(TpcFoundCluster, 6)
};

#endif
