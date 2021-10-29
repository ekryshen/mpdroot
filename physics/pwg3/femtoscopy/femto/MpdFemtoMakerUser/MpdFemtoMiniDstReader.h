/**
 * \class MpdFemtoMiniDstReader
 * \brief Reader for the MpdMiniDst format
 *
 * The reader class for MpdMiniDst format. It reads McDst and
 * converts data to the internal MpdFemtoMaker structure (MpdFemtoEvent,
 * MpdFemtoTrack, etc).
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date September 03, 2020
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoMiniDstReader_h
#define MpdFemtoMiniDstReader_h

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseEventReader.h"
// Infrastructure
#include "MpdFemtoEvent.h"
#include "MpdFemtoTrack.h"
#include "MpdFemtoV0.h"
#include "MpdFemtoKink.h"
#include "MpdFemtoXi.h"
#include "MpdFemtoEnumeration.h"
#include "MpdFemtoString.h"
#include "MpdFemtoModelHiddenInfo.h"

// MpdMiniDst headers
#include "MpdMiniDst.h"
#include "MpdMiniDstReader.h"
#include "MpdMiniEvent.h"
#include "MpdMiniTrack.h"
#include "MpdMiniMcTrack.h"
#include "MpdMiniBTofPidTraits.h"

// ROOT haders
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"

//_________________
class MpdFemtoMiniDstReader : public MpdFemtoBaseEventReader {

 public:
  /// Default constructor
  MpdFemtoMiniDstReader();
  /// Constructor
  MpdFemtoMiniDstReader(MpdMiniDstReader *miniDstReader, int debug=1);
  /// Copy constructor
  MpdFemtoMiniDstReader(const MpdFemtoMiniDstReader& copy);
  /// Assignment operator
  MpdFemtoMiniDstReader& operator=(const MpdFemtoMiniDstReader& copy);
  /// Destructor
  virtual ~MpdFemtoMiniDstReader();

  /// Return MpdFemtoEvent (back to the analysis)
  virtual MpdFemtoEvent *returnHbtEvent();
  /// Make report
  virtual MpdFemtoString report();

  /// Set miniDst reader
  void setMiniDstReader(MpdMiniDstReader *miniDstReader) { mMiniDstReader = miniDstReader; }

 private:

  /// Pointer to the MpdMiniDstReader
  MpdMiniDstReader *mMiniDstReader;
  /// Pointer to MpdFemtoEvent
  MpdFemtoEvent *mHbtEvent;

  /// Number of events passed
  long int mEventsPassed;

#ifdef __ROOT__
  ClassDef(MpdFemtoMiniDstReader, 1)
#endif
};

#endif // #define MpdFemtoMiniDstReader_h
