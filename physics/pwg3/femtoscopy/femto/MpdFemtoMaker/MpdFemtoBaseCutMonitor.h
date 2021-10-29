/**
 * \class MpdFemtoBaseCutMonitor
 * \brief The base class for cut monitors
 *
 * A cut monitor saves attributes of the entities that have passed or failed
 * the given cut. This class is the base class which is present to
 * provide a common interface for storing data.
 *
 * Cut monitors are to be used in conjunction with cut monitor handlers
 * (MpdFemtoCutMonitorHandler) - of which all standard cuts (e.g.
 * MpdFemtoEventCut) inherit from. Your cut monitor objects get added to the
 * monitorhandlers via their addCutMonitor methods, and the Fill commands get
 * called by the handler upon their fillCutMonitor method, with the particular
 * entity type.
 *
 * The default behavior of this base class is to do nothing with the
 * incoming data, as no data members are provided. It is up to the user to use
 * (or write) a subclass with relevant histograms.
 *
 * To implement a custom cut monitor, subclass this class and overload the
 * Fill method(s) corresponding to the entity type(s) you wish to monitor.
 * All other 'fill' methods should be implemented to avoid 'member hiding'
 * compiler warnings.
 *
 * All methods of this class are empty except report which returns an empty
 * string and getOutputList which returns a pointer to an empty list.
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoBaseCutMonitor_h
#define MpdFemtoBaseCutMonitor_h

// Forward declarations
class MpdFemtoEvent;
class MpdFemtoTrack;
class MpdFemtoV0;
class MpdFemtoXi;
class MpdFemtoKink;
class MpdFemtoPair;

// C++ headers
#include <iostream>

// MpdFemtoMaker headers
#include "MpdFemtoString.h"
#include "MpdFemtoParticleCollection.h"

// ROOT headers
#include "TList.h"

//_________________
class MpdFemtoBaseCutMonitor {
 public:
  /// Default constructor

  MpdFemtoBaseCutMonitor() {
    /* no-op */
  };
  /// Default destructor

  virtual ~MpdFemtoBaseCutMonitor() {
    /* no-op */
  };

  /// Report details
  virtual MpdFemtoString report();

  /// Start event processing

  virtual void eventBegin(const MpdFemtoEvent*) {
    /* no-op */
  }
  /// Finish event processing

  virtual void eventEnd(const MpdFemtoEvent*) {
    /* no-op */
  }

  /// Returns pointer to empty list
  virtual TList* getOutputList();

  /// Fill method for event
  virtual void fill(const MpdFemtoEvent*);
  /// Fill method for track
  virtual void fill(const MpdFemtoTrack*);
  /// Fill method for v0
  virtual void fill(const MpdFemtoV0*);
  /// Fill method for xi
  virtual void fill(const MpdFemtoXi*);
  /// Fill method for kink
  virtual void fill(const MpdFemtoKink*);
  /// Fill method for pair
  virtual void fill(const MpdFemtoPair*);
  /// Fill method for particle collection
  virtual void fill(const MpdFemtoParticleCollection*);
  /// Fill method for event and particle collection
  virtual void fill(const MpdFemtoEvent*, const MpdFemtoParticleCollection*);
  /// Fill method for two particle collections
  virtual void fill(const MpdFemtoParticleCollection*, const MpdFemtoParticleCollection*);

  /// Finish
  virtual void finish();
  /// Initialize cut monitor
  virtual void init();

  ClassDef(MpdFemtoBaseCutMonitor, 0)
};

//_________________
inline TList* MpdFemtoBaseCutMonitor::getOutputList() {
  TList *mOutputList = new TList();
  return mOutputList;
}

//_________________
inline MpdFemtoString MpdFemtoBaseCutMonitor::report() {
  MpdFemtoString defReport("*** no user defined fill(const MpdFemtoEvent*), take from base class");
  return defReport;
}

#endif //#define MpdFemtoBaseCutMonitor_h
