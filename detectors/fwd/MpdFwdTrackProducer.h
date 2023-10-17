//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFwdHitProducer
///
/// \brief
/// \author Evgeny Kryshen (PNPI, Gatchina)
//------------------------------------------------------------------------------------------------------------------------
#ifndef __HH_MPDFWDTRACKPRODUCER_H
#define __HH_MPDFWDTRACKPRODUCER_H

#include "FairTask.h"
class TClonesArray;

class MpdFwdTrackProducer : public FairTask {
public:
  MpdFwdTrackProducer();
  virtual ~MpdFwdTrackProducer();
  virtual InitStatus Init();
  virtual void       Exec(Option_t *opt);
  virtual void       Finish();
private:
  TClonesArray *fMcTracks;   //!
  TClonesArray *fFwdPoints;   //!
  TClonesArray *fFwdHits;     //!
  TClonesArray *fFwdTracks;   //!

  ClassDef(MpdFwdTrackProducer, 2)
};

#endif // #ifndef __HH_MPDFWDHITPRODUCER_H
