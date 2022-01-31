//
// MpdMiniMcTrack holds information about Monte Carlo particle
//

// MpdMiniDst headers
#include "MpdMiniMessMgr.h"
#include "MpdMiniMcTrack.h"

ClassImp(MpdMiniMcTrack)

   //_________________
   MpdMiniMcTrack::MpdMiniMcTrack()
   : TObject(), fId(0), fPdgId(0), fPx(0), fPy(0), fPz(0), fEnergy(0), fX(0), fY(0), fZ(0), fT(0), fMotherId(-1)
{
   // Default constructor
   if (!fRecoTrackIds.empty()) {
      fRecoTrackIds.clear();
   }
}

//_________________
MpdMiniMcTrack::MpdMiniMcTrack(const MpdMiniMcTrack &copy) : TObject()
{
   // Copy constructor
   fId    = copy.fId;
   fPdgId = copy.fPdgId;
   /* fParentIndex = copy.fParentIndex; */
   fPx           = copy.fPx;
   fPy           = copy.fPy;
   fPz           = copy.fPz;
   fEnergy       = copy.fEnergy;
   fX            = copy.fX;
   fY            = copy.fY;
   fZ            = copy.fZ;
   fT            = copy.fT;
   fRecoTrackIds = copy.fRecoTrackIds;
   fMotherId     = copy.fMotherId;
}

//_________________
MpdMiniMcTrack::~MpdMiniMcTrack()
{
   // Destructor
   /* empty */
}

//________________
void MpdMiniMcTrack::addGlobalTrackId(UShort_t id)
{

   if (!fRecoTrackIds.empty()) {
      // Assume that the new trigger is not in the list
      Bool_t isUsed = false;

      // Loop over the trigger list
      for (UInt_t iIter = 0; iIter < fRecoTrackIds.size(); iIter++) {

         // Compare triggers
         if (fRecoTrackIds.at(iIter) == id) {
            isUsed = true;
         }
      } //(UInt_t iIter=0; iIter<fRecoTrackIds.size(); iIter++)

      // If the trigger not in the list then add it
      if (!isUsed) {
         fRecoTrackIds.push_back(id);
      }
   } else {
      fRecoTrackIds.push_back(id);
   }
}

//_________________
void MpdMiniMcTrack::Print(const Char_t *option __attribute__((unused))) const
{
   LOG_INFO << "id: " << fId << " PDG code: " << fPdgId << "fourMomentum (px,py,pz,e): (" << fPx << "," << fPy << ","
            << fPz << "," << fEnergy << ")\n"
            << "fourCoordinate (x,y,z,t): (" << fX << "," << fY << "," << fZ << "," << fT << ")\n"
            << endm;
}
