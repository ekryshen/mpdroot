//------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------
// -----                     MpdStsPoint header file                   -----
// -------------------------------------------------------------------------

#ifndef MPDSTSPOINT_H
#define MPDSTSPOINT_H

#include "TObject.h"
#include "TVector3.h"
#include "FairMCPoint.h"

using namespace std;

//------------------------------------------------------------------------------------------------------------------------
class MpdStsPoint : public FairMCPoint {

public:
   /** Constructor with arguments
    *@param trackID  Index of MCTrack
    *@param detID    Detector ID
    *@param pos      Ccoordinates at entrance to active volume [cm]
    *@param mom      Momentum of track at entrance [GeV]
    *@param tof      Time since event start [ns]
    *@param length   Track length since creation [cm]
    *@param eLoss    Energy deposit [GeV]
    **/
   MpdStsPoint(Int_t trackID, Int_t detID, TVector3 posIn, TVector3 momIn, TVector3 posOut, Double_t tof,
               Double_t length, Double_t eLoss);
   MpdStsPoint(const MpdStsPoint &point) { *this = point; }

   MpdStsPoint();
   virtual ~MpdStsPoint();

   // Output to screen
   virtual void Print(const Option_t *opt) const;

   Double_t GetXout() const { return fXout; }
   Double_t GetYout() const { return fYout; }
   Double_t GetZout() const { return fZout; }
   void     PositionOut(TVector3 &pos) { pos.SetXYZ(fXout, fYout, fZout); }
   void     SetPositionOut(const TVector3 &pos);

protected:
   Double32_t fXout, fYout, fZout; // Point coordinates at exit [cm]

   ClassDef(MpdStsPoint, 1)
};

inline void MpdStsPoint::SetPositionOut(const TVector3 &pos)
{
   fXout = pos.X();
   fYout = pos.Y();
   fZout = pos.Z();
}

//------------------------------------------------------------------------------------------------------------------------
#endif
