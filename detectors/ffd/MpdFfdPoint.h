//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFfdPoint
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------

#ifndef __HH_MPDFFDPOINT_H
#define __HH_MPDFFDPOINT_H

#include <TVector3.h>
#include <TLorentzVector.h>

#include "FairMCPoint.h"
//------------------------------------------------------------------------------------------------------------------------
class MpdFfdPoint : public FairMCPoint {
public:
   enum FFDPointMode { kCherenkovPhoton = 0, kPhotoElectron };
   typedef std::vector<std::pair<double, double>> TOpContainer;

   static FFDPointMode fCurrentMode;
   /**
    *@param tid  	Index of MCTrack
    *@param duid    	Detector ID
    *@param pos      	Coordinates at entrance to active volume [cm]
    *@param mom      	Momentum of track at entrance [GeV]
    *@param tof      	Time since event start [ns]
    *@param length   	Track length since creation [cm]
    *@param eLoss    	Energy deposit [GeV]
    **/
   MpdFfdPoint(Int_t tid, Int_t duid);
   MpdFfdPoint();
   virtual ~MpdFfdPoint();

   virtual void Print(const Option_t *opt) const;

   bool AddOp(double energy, double time)
   {
      fData.push_back(std::make_pair(energy, time));
      return isClosed;
   };
   void SaveParentTrackParams(const TVector3 &posIn, const TVector3 &posOut, const TLorentzVector &P, Double_t time,
                              Double_t length);
   bool IsSame(Int_t tid, Int_t duid) const { return (GetTrackID() == tid && GetDetectorID() == duid); }
   bool IsClosed() const { return isClosed; }

   void GetPosition(TVector3 &enter, TVector3 &exit)
   {
      enter = TVector3(fX, fY, fZ);
      exit  = fPositionOut;
   }
   double              GetBeta() const { return fBeta; }
   size_t              GetEntries() const { return fData.size(); }
   const TOpContainer &GetData() const { return fData; }
   FFDPointMode        GetMode() const { return fMode; }

private:
   FFDPointMode fMode = kPhotoElectron;
   TOpContainer fData =
      {}; // pair<energy, time> [eV, ns] parameters for cherenkov photon or photo electron(dependent from fMode flag)
   TVector3 fPositionOut = {}; // [cm]
   Double_t fBeta        = 0.;
   bool     isClosed     = false;

   ClassDef(MpdFfdPoint, 2);
};
//------------------------------------------------------------------------------------------------------------------------
#endif // #ifndef __HH_MPDFFDPOINT_H
