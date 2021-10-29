/**
 * \class MpdFemtoBaseAnalysis
 * \brief The pure virtual base class for analysis
 *
 * The MpdFemtoBaseAnalysis class is a base class for analysis
 * classes must inherit from this one
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoBaseAnalysis_h
#define MpdFemtoBaseAnalysis_h

// MpdFemtoMaker headers
#include "MpdFemtoTypes.h"

// ROOT headers
#include "TList.h"
#include "TObjString.h"

// Forward declaration
class MpdFemtoEvent;

//_________________
class MpdFemtoBaseAnalysis {
 public:
  /// Default constructor

  MpdFemtoBaseAnalysis() {
    /* noop */
  }
  /// Default destructor

  virtual ~MpdFemtoBaseAnalysis() {
    /* noop */
  }

  /// Returns reports of all cuts applied and
  /// correlation functions being done
  virtual MpdFemtoString report() = 0; //!<

  /// Return list of cut settings for the analysis
  virtual TList *listSettings() = 0; //!<

  /// Obtain number of objects to be written as an output
  virtual TList *getOutputList() = 0; ///<

  /// Main machinery
  virtual void processEvent(const MpdFemtoEvent*) = 0; ///<

  /// Finish
  virtual void finish() = 0; ///<

  ClassDef(MpdFemtoBaseAnalysis, 0)
};

#endif // #define MpdFemtoBaseAnalysis_h
