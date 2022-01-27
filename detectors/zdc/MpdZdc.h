/*************************************************************************************
 *
 *         Class MpdZdc
 *
 *  Adopted for MPD by:   Elena Litvinenko
 *  e-mail:   litvin@nf.jinr.ru
 *  Version:  8-Apr-2008
 *
 ************************************************************************************/

#ifndef MPDZDC_H
#define MPDZDC_H

#include "TClonesArray.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "FairDetector.h"
#include "MpdZdcGeoPar.h"

using namespace std;

class TClonesArray;
class MpdZdcPoint;
class FairVolume;

class MpdZdc : public FairDetector {

public:
   /** Default constructor **/
   MpdZdc();

   /** Standard constructor.
    *@param name    detetcor name
    *@param active  sensitivity flag
    **/
   MpdZdc(const char *name, Bool_t active);

   /** Destructor **/
   virtual ~MpdZdc();

   /** Virtual method Initialize
    ** Initialises detector. Stores volume IDs for MUO detector and mirror.
    **/
   virtual void Initialize();

   /** Virtual method ProcessHits
    **
    ** Defines the action to be taken when a step is inside the
    ** active volume. Creates MpdZdcPoints and MpdZdcMirrorPoints and adds
    ** them to the collections.
    *@param vol  Pointer to the active volume
    **/
   virtual Bool_t ProcessHits(FairVolume *vol = 0);

   /** Virtual method EndOfEvent
    **
    ** If verbosity level is set, print hit collection at the
    ** end of the event and resets it afterwards.
    **/
   virtual void EndOfEvent();

   virtual void BeginEvent();
   /** Virtual method Register
    **
    ** Registers the hit collection in the ROOT manager.
    **/
   virtual void Register();

   /** Accessor to the hit collection **/
   virtual TClonesArray *GetCollection(Int_t iColl) const;

   /** Virtual method Print
    **
    ** Screen output of hit collection.
    **/
   virtual void Print() const;

   /** Virtual method Reset
    **
    ** Clears the hit collection
    **/
   virtual void Reset();

   /** Virtual method CopyClones
    **
    ** Copies the hit collection with a given track index offset
    *@param cl1     Origin
    *@param cl2     Target
    *@param offset  Index offset
    **/
   virtual void CopyClones(TClonesArray *cl1, TClonesArray *cl2, Int_t offset);

   /** Virtual method Construct geometry
    **
    **/
   virtual void ConstructGeometry();

   // Check whether a volume is sensitive.
   // The decision is based on the volume name. Only used in case
   // of GDML and ROOT geometry.
   // @param name    Volume name
   // @value         kTRUE if volume is sensitive, else kFALSE
   virtual Bool_t CheckIfSensitive(std::string name);

   MpdZdcPoint *GetHit(Int_t i) const;
   MpdZdcPoint *GetHit(Int_t vsc, Int_t mod, Int_t zdc) const;
   // MpdZdcPoint* GetHitH(Int_t vsch, Int_t modh) const;
   Int_t GetVSCVolId() { return fVSCVolId; }
   Int_t GetVSCHVolId() { return fVSCHVolId; }

   MpdZdcPoint *AddHit(Int_t trackID, Int_t detID, Int_t copyNo, Int_t copyNoMother, Int_t copyNoZdc, TVector3 pos,
                       TVector3 mom, Double_t tof, Double_t length, Double_t eLoss);

private:
   Int_t          fTrackID;  //!  track index
   Int_t          fVolumeID; //!  volume id
   Int_t          fEventID;  //!  event id
   TLorentzVector fPos;      //!  position
   TLorentzVector fMom;      //!  momentum
   Double32_t     fTime;     //!  time
   Double32_t     fLength;   //!  length
   Int_t          fHitNb;
   Int_t          fVSCVolId;
   Int_t          fVSCHVolId;
   Double32_t     fELoss;      //!  energy loss
   Int_t          fPosIndex;   //!
   Int_t          volDetector; //!  MC volume ID of MUO

   TClonesArray *fZdcCollection; //! Hit collection

   // reset all parameters
   void ResetParameters();

   ClassDef(MpdZdc, 2);
};

//------------------------------------------------------------------------------------------------------------------------
inline void MpdZdc::ResetParameters()
{
   fTrackID = fVolumeID = 0;
   fPos.SetXYZM(0.0, 0.0, 0.0, 0.0);
   fMom.SetXYZM(0.0, 0.0, 0.0, 0.0);
   fTime = fLength = fELoss = 0;
   fPosIndex                = 0;
};
//------------------------------------------------------------------------------------------------------------------------
#endif
