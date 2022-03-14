#ifndef MPDITSKALMANTRACK_H
#define MPDITSKALMANTRACK_H

/// \ingroup rec
/// \class MpdItsKalmanTrack
/// \brief Kalman track in MPD ITS
///
/// \author Alexander Zinchenko, LHEP JINR Dubna

#include "MpdKalmanTrack.h"
#include "MpdKalmanHit.h"
#include "MpdVector.h"
#include <TClonesArray.h>
class MpdTpcKalmanTrack;
class MpdEctKalmanTrack;
class MpdVector;
class TVector3;

class MpdItsKalmanTrack : public MpdKalmanTrack {

public:
   MpdItsKalmanTrack();                                                                         ///< Default ctor
   virtual ~MpdItsKalmanTrack();                                                                ///< Destructor
   MpdItsKalmanTrack(MpdKalmanHit *hitOut, MpdKalmanHit *hitIn, TVector3 &vertex, Double_t pt); ///< Ctor from 2 hits
   MpdItsKalmanTrack(const MpdItsKalmanTrack &track);                                           ///< copy constructor
   MpdItsKalmanTrack &operator=(const MpdItsKalmanTrack &track);                                // assignment operator
   MpdItsKalmanTrack(const MpdTpcKalmanTrack &track); ///< constructor from TPC track

   MpdItsKalmanTrack(const MpdVector &track1,
                     TVector3        &vec); ///< constructor from CellTrack track /// MpdVector&MpdCellTrack&

   Int_t         GetNofTrHits() const { return fTrHits->GetEntriesFast(); } ///< get number of track hits
   TClonesArray *GetTrHits() const { return fTrHits; }                      ///< get track hits

   // Int_t GetNofKHits() const { return fKHits->GetEntriesFast(); } //attention!!!
   // TClonesArray *GetKHits() const { return fKHits; } // attention!!!

   Int_t    GetNofIts() const { return fNofIts; }   ///< get number of ITS hits
   Double_t GetChi2Its() const { return fChi2Its; } ///< get Chi2 of ITS fit

   void StartBack(); ///< prepare for back tracing

   void SetNofIts(Int_t nhits) { fNofIts = nhits; }    ///< set number of ITS hits
   void SetChi2Its(Double_t chi2) { fChi2Its = chi2; } ///< set Chi2 of ITS fit

   Bool_t IsSortable() const { return kTRUE; }
   Int_t  Compare(const TObject *track) const; ///< sort in descending order in Pt
   void   Reset();                             ///< reset track (similar to destructor)

private:
   void EvalCovar(const MpdKalmanHit *hitOut, const MpdKalmanHit *hitIn, Double_t *posOut,
                  Double_t *posIn); ///< evaluate covar. matrix

   //  TClonesArray *fKHits; //attention!!!
   TClonesArray *fTrHits;  ///< track hits
   Int_t         fNofIts;  ///< number of ITS hits
   Double32_t    fChi2Its; ///< Chi2 of the fit in ITS

   ClassDef(MpdItsKalmanTrack, 1);
};
#endif
