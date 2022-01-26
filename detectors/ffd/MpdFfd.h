//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFfd
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------

#ifndef __HH_MPDFFD_H
#define __HH_MPDFFD_H

#include <iostream>

#include <vector>
#include <map>

#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "FairDetector.h"

#include "MpdFfdPoint.h"

// class MpdFfdPoint;
class FairVolume;

//------------------------------------------------------------------------------------------------------------------------
class MpdFfd : public FairDetector {
public:
   // *@param name    detector name
   // *@param active  sensitivity flag
   MpdFfd(const char *name, Bool_t active, Int_t verbose = 1);
   MpdFfd();
   virtual ~MpdFfd();

   // Defines the action to be taken when a step is inside the
   // active volume. Creates MpdFfdPoints and adds them to the collection.
   // @param vol  Pointer to the active volume
   virtual Bool_t ProcessHits(FairVolume *vol = 0);

   // If verbosity level is set, print hit collection at the
   // end of the event and resets it afterwards.
   virtual void EndOfEvent();

   // Registers the hit collection in the ROOT manager.
   virtual void Register();

   // Accessor to the hit collection
   virtual TClonesArray *GetCollection(Int_t iColl) const;

   // Screen output of hit collection.
   virtual void Print() const;

   // Clears the hit collection
   virtual void Reset();

   // *@param cl1     Origin
   // *@param cl2     Target
   // *@param offset  Index offset
   virtual void CopyClones(TClonesArray *cl1, TClonesArray *cl2, Int_t offset);

   // Constructs the FFD geometry
   virtual void ConstructGeometry();

   // Construct the geometry from an ASCII geometry file
   virtual void ConstructAsciiGeometry();

   // Check whether a volume is sensitive.
   // The decision is based on the volume name. Only used in case
   // of GDML and ROOT geometry.
   // @param name    Volume name
   // @value         kTRUE if volume is sensitive, else kFALSE
   virtual Bool_t CheckIfSensitive(std::string name);

private:
   TLorentzVector fPos = {}; //!  position
   TLorentzVector fMom = {}; //!  momentum

   // cached looking for
   MpdFfdPoint *cPoint = nullptr;
   Int_t        cPtid = -1, cSuid = -1;

   Int_t         fPosIndex;  //!
   TClonesArray *aFfdPoints; //! Hit collection

   MpdFfdPoint *FindPoint(Int_t ptid, Int_t suid);
   MpdFfdPoint *CreatePoint(Int_t ptid, Int_t suid);

   // Resets the private members for the track parameters
   void ResetParameters();

   struct TrackParam {
      TrackParam(Int_t id, const TVector3 &in, const TLorentzVector &P, Double_t t,
                 Double_t l) // init params at track enter to SV
      {
         suid   = id;
         posIn  = in;
         mom    = P;
         time   = t;
         length = l;
      }

      Int_t          suid  = -1;
      TVector3       posIn = {}, posOut = {};
      TLorentzVector mom  = {};
      Double_t       time = 0., length = 0.;
   };

   typedef std::multimap<Int_t, TrackParam> Tparams;

   Tparams params; //!   key = tid

   std::pair<Tparams::iterator, bool> FindTrackParams(Int_t tid, Int_t suid);

   ClassDef(MpdFfd, 3)
};
//------------------------------------------------------------------------------------------------------------------------
inline void MpdFfd::ResetParameters()
{
   fPos = fMom = {};
   fPosIndex   = {};
};
//------------------------------------------------------------------------------------------------------------------------
#endif // #ifndef __HH_MPDFFD_H
