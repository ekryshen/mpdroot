//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFwd
///
/// \brief
/// \author Evgeny Kryshen (PNPI, Gatchina)
//------------------------------------------------------------------------------------------------------------------------
#ifndef __HH_MPDFWD_H
#define __HH_MPDFWD_H

#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "FairDetector.h"
class FairVolume;

//------------------------------------------------------------------------------------------------------------------------
class MpdFwd : public FairDetector {
public:
   MpdFwd(const char *name="FWD");
   virtual ~MpdFwd();
   // Defines the action to be taken when a step is inside the
   // active volume. Creates MpdFwdPoints and adds them to the collection.
   // @param vol  Pointer to the active volume
   virtual Bool_t ProcessHits(FairVolume *vol = 0);
   // If verbosity level is set, print hit collection at the
   // end of the event and resets it afterwards.
   virtual void EndOfEvent();
   // Registers the hit collection in the ROOT manager.
   virtual void Register();
   // Accessor to the hit collection
   virtual TClonesArray *GetCollection(Int_t iColl) const;
   // Clears the hit collection
   virtual void Reset();
   // Constructs the FWD geometry
   virtual void ConstructGeometry();
private:
  TClonesArray *fFwdPoints; //! McPoint collection
  Int_t fTrackID;           //!  track index
  TLorentzVector fPosIn;    //!  position in
  TLorentzVector fPosOut;   //!  position out
  Double32_t fTime;         //!  time
  Double32_t fLength;       //!  length
  Double32_t fELoss;        //!  energy loss

  ClassDef(MpdFwd, 1)
};

#endif // #ifndef __HH_MPDFWD_H
