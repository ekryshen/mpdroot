#ifndef MPD_ITSTOTPCMATCHING_H
#define MPD_ITSTOTPCMATCHING_H

/// \ingroup rec
/// \class MpdItsToTpcMatching
/// \brief Track finder in MPD Inner Tracking System (ITS) using cellular automaton (CA)
///
/// \author Dmitrii Zinchenko, Alexander Zinchenko, LHEP JINR Dubna

#include "MpdItsHit5spd.h"
#include "MpdVector.h"

#include "FairTask.h"

#include "TH1.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TVector2.h"

#include <map>
#include <set>

class MpdKalmanHit;
class MpdKalmanTrack;
class MpdItsKalmanTrack;
class MpdTpcKalmanTrack;

class MpdItsToTpcMatching : public FairTask {
public:
   /** Constructor **/
   MpdItsToTpcMatching(const char *name = "MpdVectorFinder", Int_t iVerbose = 1);

   /** Destructor **/
   virtual ~MpdItsToTpcMatching();
   /// * FairTask methods

   /** Intialisation at begin of run. To be implemented in the derived class.
    *@value  Success   If not kSUCCESS, task will be set inactive.
    **/
   InitStatus Init();

   /** Reinitialisation.
    *@value  Success   If not kSUCCESS, task will be set inactive.
    **/
   InitStatus ReInit();

   /** Intialise parameter containers.
    **/
   void SetParContainers();

   void GetItsTracks(std::multimap<Float_t, MpdItsKalmanTrack *> &m, std::multimap<Float_t, MpdItsKalmanTrack *> &mPhi);
   void RefitItsTo27(std::multimap<Float_t, MpdItsKalmanTrack *> &m, std::multimap<Float_t, MpdItsKalmanTrack *> &mPhi);
   void RefitTpcTo27(std::multimap<Float_t, MpdTpcKalmanTrack *> &m, std::multimap<Float_t, MpdTpcKalmanTrack *> &mPhi);
   void Exec(Option_t *option);

   /** Action after each event. **/
   void Finish();
   void Reset();
   // void FillGeoScheme(); // fill Kalman filter geometry manager   ///???

private:
   // Int_t fNPass;
   Int_t fExact; //!< exact ID match if != 0

   Int_t           fNTracks;        // number of found tracks 10.02
   TClonesArray   *fItsPoints;      //! ITS MC points
   TClonesArray   *fItsTracks;      //! ITS tracks
   TClonesArray   *fTpcTracks;      //!< array of TPC tracks
   TClonesArray   *fMCTracks;       //!< array of MC tracks
   TClonesArray   *fTracks;         //! array of ITS tracks
   TClonesArray   *fTracks1;        //! array of TPC tracks
   TClonesArray   *fTracksRefit;    //! array of refit ITS tracks
   TClonesArray   *fTpcTracksRefit; //! array of refit TPC tracks;
   TClonesArray   *fKHits;          //! array of Kalman hits from Vector Finder
   std::set<Int_t> fTpcIndSet;      //! unmatched TPC track indices

   ClassDef(MpdItsToTpcMatching, 0);
};

#endif
