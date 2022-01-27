#ifndef MPDTPCKALMANTRACK_H
#define MPDTPCKALMANTRACK_H

/// \ingroup rec
/// \class MpdTpcKalmanTrack
/// \brief Kalman track in MPD TPC
///
/// \author Alexander Zinchenko, LHEP JINR Dubna

#include "MpdKalmanTrack.h"
#include "MpdKalmanHit.h"
#include <TClonesArray.h>
class TVector3;

class MpdTpcKalmanTrack : public MpdKalmanTrack {

public:
   MpdTpcKalmanTrack();                               ///< Default ctor
   MpdTpcKalmanTrack(Double_t pos, TVector3 &vertex); // ctor used by KFV0 finder
   virtual ~MpdTpcKalmanTrack();                      ///< Destructor
   MpdTpcKalmanTrack(MpdKalmanHit *hitOut, MpdKalmanHit *hitIn, TVector3 &vertex, Double_t pt,
                     Double_t *params); ///< Ctor from 2 hits
   MpdTpcKalmanTrack(MpdKalmanHit *hitOut, MpdKalmanHit *hitIn, TVector3 &vertex, TVector3 &posOut, TVector3 &posIn,
                     Double_t pt, Double_t *params);                             ///< Ctor from 2 hits and 2 TVector3
   MpdTpcKalmanTrack(const MpdTpcKalmanTrack &track);                            ///< copy constructor
   MpdTpcKalmanTrack &operator=(const MpdTpcKalmanTrack &track);                 // assignment operator
   Int_t              GetNofTrHits() const { return fTrHits->GetEntriesFast(); } ///< get number of track hits
   TClonesArray *     GetTrHits() const { return fTrHits; }                      ///< get track hits

   void StartBack(); ///< prepare for back tracing

   Bool_t IsSortable() const { return kTRUE; }
   Int_t  Compare(const TObject *track) const; ///< sort in descending order in Pt
   void   Reset();                             ///< reset track (similar to destructor)
   Bool_t GetRecoQuality(Double_t dist       = 1.5,
                         Double_t percentage = 0.5); ///< returns kTRUE if number of hits closer to boundaries than dist
                                                     ///< divided by nHits is larger than percentage

private:
   void EvalParams(MpdKalmanHit *hitOut, const MpdKalmanHit *hitIn, Double_t *parOut, Double_t *parIn, Double_t pt,
                   Double_t *params); ///< evaluate track params
   void EvalCovar(const MpdKalmanHit *hitOut, const MpdKalmanHit *hitIn, Double_t *parOut,
                  Double_t *parIn); ///< evaluate covar. matrix

   TClonesArray *fTrHits; ///< track hits

   ClassDef(MpdTpcKalmanTrack, 1);
};
#endif
