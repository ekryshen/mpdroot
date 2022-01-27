//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_ETOF_H
#define __MPD_ETOF_H 1

#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "MpdTof.h"
#include "MpdTofUtils.h"

#include "FairDetector.h"

class MpdTofPoint;
class FairVolume;
//------------------------------------------------------------------------------------------------------------------------
class MpdEtof : public FairDetector {
   // Track information to be stored until the track leaves the active volume.
   Int_t          fTrackID;  //!  track index
   Int_t          fVolumeID; //!  volume id
   TLorentzVector fPos;      //!  position
   TLorentzVector fMom;      //!  momentum
   Double_t       fTime;     //!  time
   Double_t       fLength;   //!  length
   Double_t       fELoss;    //!  energy loss

   Int_t         fPosIndex; //!
   TClonesArray *aTofHits;  //! Hit collection
   const double  nan;

   MpdTofPoint *AddPoint(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom, Double_t time, Double_t length,
                         Double_t eLoss);
   void         ResetParameters();

public:
   MpdEtof(const char *name = "ETOF", Bool_t active = kTRUE);
   virtual ~MpdEtof();

   virtual Bool_t        ProcessHits(FairVolume *vol = 0);
   virtual void          EndOfEvent();
   virtual void          Register();
   virtual TClonesArray *GetCollection(Int_t iColl) const;
   virtual void          Print() const;
   virtual void          Reset();
   virtual void          CopyClones(TClonesArray *cl1, TClonesArray *cl2, Int_t offset);
   virtual void          ConstructGeometry();

   ClassDef(MpdEtof, 3);
};
//------------------------------------------------------------------------------------------------------------------------
inline void MpdEtof::ResetParameters()
{
   fTrackID = fVolumeID = 0;
   fPos.SetXYZM(nan, nan, nan, nan);
   fMom.SetXYZM(nan, nan, nan, nan);
   fTime = fLength = fELoss = nan;
   fPosIndex                = 0;
};
//------------------------------------------------------------------------------------------------------------------------
#endif
