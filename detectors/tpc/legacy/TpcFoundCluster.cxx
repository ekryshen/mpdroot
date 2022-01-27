//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class TpcFoundCluster
//
//
// Environment:
//      Software developed for the MPD at NICA.
//
// Author List:
//      Salmin Roman           (original author)
//
//
//-----------------------------------------------------------

#include "TpcFoundCluster.h"
#include "TpcPadPlane.h"

ClassImp(TpcFoundCluster);

TpcFoundCluster::TpcFoundCluster() : _localPos(0, 0, 0), _sectorID(-1) {}

TpcFoundCluster::TpcFoundCluster(const TVector3 &pos, Int_t sectorID, UInt_t peak, Double_t valey, Double_t charge,
                                 UInt_t MCTrackID)
   : _localPos(pos), _sectorID(sectorID), _peak(peak), _meanValey(valey), _charge(charge),
     _MCTrackID(MCTrackID) // charge added by Sergey Merts
{
}

TpcFoundCluster::~TpcFoundCluster() {}

TVector3 TpcFoundCluster::globalPos() const
{
   const TpcPadPlane *padPlane = TpcPadPlane::Instance();
   return padPlane->fromSectorReferenceFrame(_localPos, _sectorID);
}
