//--------------------------------------------------------------------
// Description:
//      Test-driven interface (Abstract base class) for TPC Hits
//
//      TPC Hit is defined by location within its cluster
//
// Environment:
//      Software for the MPD Detector at NICA.
//
// Authors:
//      Slavomir Hnatic
//      JINR, December, 2022
//--------------------------------------------------------------------

#ifndef ABSTRACTTPCHIT_HH
#define ABSTRACTTPCHIT_HH

#include <TVector3.h>
#include <FairHit.h>

class AbstractTpcHit : public FairHit {

public:
   AbstractTpcHit() {}
   // Constructor for compatibility with MpdTpcHit. Hit parameters are passed on to Fairhit constructor.
   AbstractTpcHit(int detID, TVector3 pos, TVector3 dpos, int index) : FairHit(detID, pos, dpos, index) {}
   virtual ~AbstractTpcHit() {}

   virtual int    GetClusterID() const         = 0;
   virtual double GetPadCoordinate() const     = 0;
   virtual double GetTimeBinCoordinate() const = 0;
   virtual float  GetTotalSignal() const       = 0;

private:
   ClassDef(AbstractTpcHit, 1);
};

#endif
