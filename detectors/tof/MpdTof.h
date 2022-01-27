//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_TOF_H
#define __MPD_TOF_H

//------------------------------------------------------------------------------------------------------------------------
/// \class MpdTof, version = 8 (14 sectors, 20 detectors, 72 strips)
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------

#include <TLorentzVector.h>
#include <TVector3.h>

#include "MpdTofUtils.h"

#include "FairDetector.h"
//------------------------------------------------------------------------------------------------------------------------
class MpdTofPoint;
class FairVolume;
class TClonesArray;
//------------------------------------------------------------------------------------------------------------------------
class MpdTof : public FairDetector {
   // Track information to be stored until the track leaves the active volume.
   Int_t          fTrackID;  //!  track index
   Int_t          fVolumeID; //!  volume id
   TLorentzVector fPos;      //!  position
   TLorentzVector fMom;      //!  momentum
   Double_t       fTime;     //!  time
   Double_t       fLength;   //!  length
   Double_t       fELoss;    //!  energy loss

   Int_t         fPosIndex = 0;       //!
   TClonesArray *aTofHits  = nullptr; //! Hit collection

   void         ConstructAsciiGeometry();
   MpdTofPoint *AddPoint(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom, Double_t time, Double_t length,
                         Double_t eLoss);
   void         ResetParameters();

public:
   MpdTof(const char *name = "TOF", Bool_t active = kTRUE);
   virtual ~MpdTof();

   virtual Bool_t        ProcessHits(FairVolume *vol = nullptr);
   virtual void          EndOfEvent();
   virtual void          Register();
   virtual TClonesArray *GetCollection(Int_t iColl) const;
   virtual void          Print() const;
   virtual void          Reset();
   virtual void          CopyClones(TClonesArray *from, TClonesArray *to, Int_t offset);
   virtual void          ConstructGeometry();
   virtual Bool_t        CheckIfSensitive(std::string name);

   static void Dump(TClonesArray *aHits, TClonesArray *aPoints, TClonesArray *aTracks, const char *comment = nullptr,
                    std::ostream &os = std::cout);
   static void Print(const TVector3 &, const char *comment = nullptr, std::ostream &os = std::cout);
   static void GetDelta(const TVector3 &mcPos, const TVector3 &estPos, double &dev, double &devZ, double &devR,
                        double &devPhi);

   ClassDef(MpdTof, 3);
};
//------------------------------------------------------------------------------------------------------------------------
#endif // #ifndef __MPD_TOF_H
