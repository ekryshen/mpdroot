// ------------------------------------------------------------------------
// -----                     MpdZdcPoint header file                  -----
// -----                     litvin@nf.jinr.ru                        -----
// -----                     Last updated 22-Feb-2012                 -----
// ------------------------------------------------------------------------

#ifndef MPDZDCPOINT_H
#define MPDZDCPOINT_H

#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "FairMCPoint.h"
#include <map>

using namespace std;

class MpdZdcPoint : public FairMCPoint {

public:
   /** Default constructor **/
   MpdZdcPoint();

   /** Constructor with arguments
    *@param trackID  Index of MCTrack
    *@param detID    Detector ID (at present, volume MC number)
    *@param copyNo         Number of active layer inside ZDC module
    *@param copyNoMother   ZDC module number
    *@param copyNoZdc
    *@param pos      Coordinates  [cm]
    *@param mom      Momentum of track [GeV]
    *@param tof      Time since event start [ns]
    *@param length   Track length since creation [cm]
    *@param eLoss    Energy deposit [GeV]
    *@param EventId
    **/

   MpdZdcPoint(Int_t trackID, Int_t detID, Int_t copyNo, Int_t copyNoMother, Int_t copyNoZdc,
               // Int_t copyNo_h, Int_t copyNoMother_h,
               TVector3 pos, TVector3 mom, Double_t tof, Double_t length, Double_t eLoss, UInt_t EventId = 0);

   /** Copy constructor **/
   MpdZdcPoint(const MpdZdcPoint &point) { *this = point; };

   /** Destructor **/
   virtual ~MpdZdcPoint();

   /** Accessors **/
   Short_t GetCopy() const { return nCopy; };
   Short_t GetCopyMother() const { return nCopyMother; }; // module
   Short_t GetCopyZdc() const { return nCopyZdc; };       // zdc (left, right)

   // Short_t GetCopyH()        const {return nCopy_h; };
   // Short_t GetCopyMotherH()  const {return nCopyMother_h; };//module with hole

   /** Modifiers **/
   void SetCopy(Short_t i) { nCopy = i; };
   void SetCopyMother(Short_t i) { nCopyMother = i; };
   void SetCopyZdc(Short_t i) { nCopyZdc = i; };

   // void SetCopyH(Short_t i)          { nCopy_h    = i; };
   // void SetCopyMotherH(Short_t i)    { nCopyMother_h  = i; };

   /** Output to screen **/
   virtual void Print(const Option_t *opt) const;

   //  fTrackID    = trackID;
   //  fDetectorID = detID;

   //  fX          = pos.X();
   //  fY          = pos.Y();
   //  fZ          = pos.Z();
   //  fPx         = mom.Px();
   //  fPy         = mom.Py();
   //  fPz         = mom.Pz();
   //  fTime       = tof;
   //  fLength     = length;
   //  fELoss      = eLoss;
   //  fEventId

   void AddVSC(Int_t trackID, Int_t detID, Int_t idvsc, Int_t idmod, Int_t idzdc, TVector3 pos, TVector3 mom,
               Double_t dt, Double_t dl, Double_t de)
   {
      if (nCopy != idvsc) cerr << "Warning: idvsc not equal in NA61VetoHit::AddVSC";
      if (nCopyMother != idmod) cerr << "Warning: idmod not equal in NA61VetoHit::AddVSC";
      if (nCopyZdc != idzdc) cerr << "Warning: idzdc not equal in NA61VetoHit::AddVSC";
      fTrackID    = trackID;
      fDetectorID = detID;
      fX          = pos.X();
      fY          = pos.Y();
      fZ          = pos.Z();
      fPx         = mom.Px();
      fPy         = mom.Py();
      fPz         = mom.Pz();
      nCopy       = idvsc;
      nCopyMother = idmod;
      nCopyZdc    = idzdc;
      fELoss += de;
      fLength = dl;
      fTime   = dt;
   }
   /*
     void AddVSCH(Int_t trackID, Int_t detID, Int_t idvsc, Int_t idmod,TVector3 pos,TVector3 mom,Double_t dt, Double_t
    dl, Double_t de) { if(nCopy_h != idvsc) cerr << "Warning: idvsc not equal in NA61VetoHit::AddVSCH"; if(nCopyMother_h
    != idmod) cerr << "Warning: idmod not equal in NA61VetoHit::AddVSCH"; fTrackID=trackID; fDetectorID=detID;
       fX=pos.X(); fY=pos.Y(); fZ=pos.Z();
       fPx=mom.Px(); fPy=mom.Py(); fPz=mom.Pz();
       nCopy_h=idvsc; nCopyMother_h=idmod; fELoss += de; fLength += dl; fTime+=dt;
     }
   */
   void IncreaseTimesVR(Float_t e, Float_t time);

   std::map<Float_t, Float_t> GetContribTimes() { return fContribTimes; }

protected:
   Short_t nCopy;       // Copy number
   Short_t nCopyMother; // Copy number of mother volume
   Short_t nCopyZdc;    // Copy number of zdc (left, right)

   // Short_t nCopy_h;                // Copy number
   // Short_t nCopyMother_h;          // Copy number of mother volume

   std::map<Float_t, Float_t> fContribTimes; // corrected time vs. deposited energy

   ClassDef(MpdZdcPoint, 4);
};
#endif
