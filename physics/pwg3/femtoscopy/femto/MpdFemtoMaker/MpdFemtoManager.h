/**
 * \class MpdFemtoManager
 * \brief The Manager is the top-level object containing an EventReader
 * (the input), and collections of EventWriters and MpdFemtoAnalyses
 * (the outputs)
 *
 * Manager is the top-level object containing an EventReader
 * (the input), and collections of EventWriters and MpdFemtoAnalyses
 * (the outputs).
 *
 * A manager object is owned by an MpdFemtoTaskFemto object which
 * calls the `processEvent()` method, which reads an MpdFemtoEvent
 * from the input files, and forwards it to the `processEvent` method
 * in each output analysis, which is responsible for carrying out the
 * actual cuts & computation.
 *
 * MpdFemtoManager objects "own" the EventReader, Analyses, and
 * EventWriters added to them, and is responsible for deleting them
 * upon its own destruction.
 *
 * MpdFemtoManager objects are not copyable, as the MpdFemtoAnalysis
 * objects they contain have no means of copying/cloning.
 * Denying copyability by making the copy constructor and assignment
 * operator private prevents potential dangling pointer (segfault)
 * errors.
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoManager_h
#define MpdFemtoManager_h

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseAnalysis.h"
#include "MpdFemtoBaseEventReader.h"
#include "MpdFemtoBaseEventWriter.h"
// Infrastructure
#include "MpdFemtoTypes.h"
#include "MpdFemtoAnalysisCollection.h"
#include "MpdFemtoEventWriterCollection.h"
#include "MpdFemtoEvent.h"

//_________________
class MpdFemtoManager {
 public:
  /// Default constructor
  MpdFemtoManager();
  /// Copy constructor
  MpdFemtoManager(const MpdFemtoManager& copy);
  /// Copy constructor
  MpdFemtoManager& operator=(const MpdFemtoManager& man);
  /// Default destructor
  virtual ~MpdFemtoManager();

  /// Return pointer to the analysis collection
  MpdFemtoAnalysisCollection *analysisCollection() {
    return mAnalysisCollection;
  }
  /// Access to the n-th analysis within Collection
  MpdFemtoBaseAnalysis *analysis(int n);
  /// Add analysis
  void addAnalysis(MpdFemtoBaseAnalysis *analysis) {
    mAnalysisCollection->push_back(analysis);
  }

  /// Return pointer to the Collection of event writers
  MpdFemtoEventWriterCollection* eventWriterCollection() {
    return mEventWriterCollection;
  }
  /// Access to n-th EventWriter within Collection
  MpdFemtoBaseEventWriter *eventWriter(int n);
  /// Add event writer
  void setEventWriter(MpdFemtoBaseEventWriter *writer) {
    addEventWriter(writer);
  }
  /// Add event writer
  void addEventWriter(MpdFemtoBaseEventWriter* writer) {
    mEventWriterCollection->push_back(writer);
  }

  /// Return event reader
  MpdFemtoBaseEventReader* eventReader() {
    return mEventReader;
  }
  /// Add EventReader
  void setEventReader(MpdFemtoBaseEventReader* reader) {
    mEventReader = reader;
  }

  /// Calls `init()` on all owned EventWriters
  /// Returns 0 for success, 1 for failure.
  int init();
  /// A "0" return value means success - otherwise quit
  int processEvent();
  /// Calls `Finish()` on the EventReader, EventWriters, and the Analyses.
  void finish();

  /// Construct report
  MpdFemtoString report();

 private:

  /// Pointer to a collection of analyses
  MpdFemtoAnalysisCollection *mAnalysisCollection;
  /// Pointer to event reader
  MpdFemtoBaseEventReader *mEventReader;
  /// Pointer to a collection of event writers
  MpdFemtoEventWriterCollection *mEventWriterCollection;

  ClassDef(MpdFemtoManager, 0)
};

#endif // MpdFemtoManager_h
