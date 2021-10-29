/**
 * \class MpdFemtoBaseEventReader
 * \brief The pure virtual base class for femto event readers
 *
 * All event readers must inherit from this one
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoBaseEventReader_h
#define MpdFemtoBaseEventReader_h

// C++ headers
#include <iostream>

// Forward declarations
class MpdFemtoEvent;
class MpdFemtoBaseEventCut;
class MpdFemtoBaseTrackCut;
class MpdFemtoBaseV0Cut;
class MpdFemtoBaseXiCut;
class MpdFemtoBaseKinkCut;

// MpdFemtoMaker headers
#include "MpdFemtoString.h"

//_________________
class MpdFemtoBaseEventReader {
 public:
  /// Default constructor
  ///
  /// Even though it's only a base class and never constructed, if you don't
  /// have an implementation, you get "AliFemtoEventReader type_info node" upon
  /// dynamical loading
  ///
  /// All pointers are set to NULL, the status is set to 0 (good), and debug is
  /// set to 1 (print debug information in methods which run once)
  ///
  MpdFemtoBaseEventReader();

  /// Copy constructor
  ///
  /// This performs a shallow copy, so both the origial and new event readers
  /// point to the same cut objects.
  MpdFemtoBaseEventReader(const MpdFemtoBaseEventReader& copy);

  /// Assignment Operator
  /// Performs shallow copy of members
  MpdFemtoBaseEventReader& operator=(const MpdFemtoBaseEventReader& copy);

  /// Destructor
  ///
  /// No members are deleted - it is up to the entity creating the cuts to
  /// delete them after the event reader has run its course
  virtual ~MpdFemtoBaseEventReader() {
    /* empty */
  }

  /// Concrete subclasses MUST implement this method, which creates the MpdFemtoEvent
  virtual MpdFemtoEvent* returnHbtEvent() = 0;

  /// User-written method to return string describing reader
  /// including whatever "early" cuts are being done
  virtual MpdFemtoString report();

  /// Next method does NOT need to be implemented, in which case the
  /// "default" method below is executed
  virtual int writeHbtEvent(MpdFemtoEvent* /* event */) {
    std::cout << "No WriteHbtEvent implemented" << std::endl;
    return 0;
  }

  // Next two are optional but would make sense for, e.g., opening and closing a file

  /// Initialization
  virtual int init(const char* ReadWrite, MpdFemtoString& Message);

  /// Finalization
  virtual void finish() {
    /* empty */
  }

  /// MpdFemtoManager looks at this for guidance if it gets null pointer from ReturnHbtEvent
  int status() {
    return mReaderStatus;
  }

  /// Set event cut
  virtual void setEventCut(MpdFemtoBaseEventCut* ecut) {
    mEventCut = ecut;
  }
  
  /// Set track cut
  virtual void setTrackCut(MpdFemtoBaseTrackCut* pcut) {
    mTrackCut = pcut;
  }
  
  /// Set v0 cut
  virtual void setV0Cut(MpdFemtoBaseV0Cut* pcut) {
    mV0Cut = pcut;
  }
  
  /// Set xi cut
  virtual void setXiCut(MpdFemtoBaseXiCut* pcut) {
    mXiCut = pcut;
  }
  
  /// Set kink cut
  virtual void setKinkCut(MpdFemtoBaseKinkCut* pcut) {
    mKinkCut = pcut;
  }

  /// Return pointer to event cut
  virtual MpdFemtoBaseEventCut* eventCut() {
    return mEventCut;
  }
  
  /// Return pointer to track cut
  virtual MpdFemtoBaseTrackCut* trackCut() {
    return mTrackCut;
  }
  
  /// Return pointer to V0 cut
  virtual MpdFemtoBaseV0Cut* v0Cut() {
    return mV0Cut;
  }
  
  /// Return pointer to V0 cut
  virtual MpdFemtoBaseXiCut* xiCut() {
    return mXiCut;
  }
  
  /// Return pointer to kink cut
  virtual MpdFemtoBaseKinkCut* kinkCut() {
    return mKinkCut;
  }

  /**
   * Controls the amount of debug information printed.
   * The code indicates which functions should print debug statements:
   *
   * 0: no output at all
   * 1: once (e.g. in constructor, finsh
   * 2: once per event
   * 3: once per track
   * 4: once per pair
   */
  int debug() {
    return mDebug;
  }

  /**
   * Set debug level:
   * 0: no output at all
   * 1: once (e.g. in constructor, finsh
   * 2: once per event
   * 3: once per track
   * 4: once per pair
   */
  void setDebug(int d) {
    mDebug = d;
  }

 protected:
  /// Link to the front-loaded event cut
  MpdFemtoBaseEventCut *mEventCut; //!<!
  /// Link to the front-loaded track cut
  MpdFemtoBaseTrackCut *mTrackCut; //!<!
  /// Link to the front-loaded V0 cut
  MpdFemtoBaseV0Cut *mV0Cut; //!<!
  /// Link to the front-loaded Xi cut
  MpdFemtoBaseXiCut *mXiCut; //!<!
  /// Link to the front-loaded Kink cut
  MpdFemtoBaseKinkCut *mKinkCut; //!<!
  /// status: 0 - good
  int mReaderStatus; ///<
  /**
   * Debug level:
   * 0: no output at all
   * 1: once (e.g. in constructor, finsh
   * 2: once per event
   * 3: once per track
   * 4: once per pair
   */
  int mDebug; ///<

  ClassDef(MpdFemtoBaseEventReader, 0)
};

#endif // #define MpdFemtoBaseEventReader_h
