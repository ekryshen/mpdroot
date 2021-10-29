/**
 * \class MpdFemtoMultiTrackCut
 * \brief Allows to add several cuts to the list
 *
 * The MpdFemtoMultiTrackCut class allows to add many track cuts
 * into the list
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoMultiTrackCut_h
#define MpdFemtoMultiTrackCut_h

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseTrackCut.h"
// Infrastructure
#include "MpdFemtoTypes.h"
#include "MpdFemtoTrack.h"
#include "MpdFemtoV0.h"
#include "MpdFemtoTrackCutCollection.h"

//_________________
class MpdFemtoMultiTrackCut : public MpdFemtoBaseTrackCut {
 public:
  /// Default constructor
  MpdFemtoMultiTrackCut();
  /// Copy constructor
  MpdFemtoMultiTrackCut(const MpdFemtoMultiTrackCut& copy);
  /// Assignment operator
  MpdFemtoMultiTrackCut& operator=(const MpdFemtoMultiTrackCut& copy);
  /// Destructor
  virtual ~MpdFemtoMultiTrackCut();

  /// User-written method to return string describing cuts
  virtual MpdFemtoString report();
  /// True if passes, false if not
  virtual bool pass(const MpdFemtoTrack* track);

  /// Add track cut
  virtual void addTrackCut(MpdFemtoBaseTrackCut*);
  /// Start event
  virtual void eventBegin(const MpdFemtoEvent*);
  /// Finish event
  virtual void eventEnd(const MpdFemtoEvent*);

  /// Return track type
  MpdFemtoParticleType type() {
    return hbtTrack;
  }

  /// Clone cut
  virtual MpdFemtoMultiTrackCut* clone();

 private:
  /// Pointer to a track cut collection
  MpdFemtoTrackCutCollection* mCutCollection;

  ClassDef(MpdFemtoMultiTrackCut, 0)
};

#endif // #define MpdFemtoMultiTrackCut_h
