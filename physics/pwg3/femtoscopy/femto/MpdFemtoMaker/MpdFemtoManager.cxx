//
// MpdFemtoManager: main class managing femtoscopic analysis
//

// C/C++ headers
#include <cstdio>

// MpdFemtoMaker headers
#include "MpdFemtoManager.h"

ClassImp(MpdFemtoManager)

//_________________
MpdFemtoManager::MpdFemtoManager() :
mAnalysisCollection(nullptr),
mEventReader(nullptr),
mEventWriterCollection(nullptr) {
  // Constructor
  mAnalysisCollection = new MpdFemtoAnalysisCollection;
  mEventWriterCollection = new MpdFemtoEventWriterCollection;
}

//_________________
MpdFemtoManager::MpdFemtoManager(const MpdFemtoManager &copy) :
mAnalysisCollection(new MpdFemtoAnalysisCollection),
mEventReader(copy.mEventReader),
mEventWriterCollection(new MpdFemtoEventWriterCollection) {
  // Copy constructor
  MpdFemtoAnalysisIterator AnalysisIter;
  for (AnalysisIter = copy.mAnalysisCollection->begin();
       AnalysisIter != copy.mAnalysisCollection->end();
       AnalysisIter++) {
    mAnalysisCollection->push_back(*AnalysisIter);
  }

  MpdFemtoEventWriterIterator EventWriterIter;
  for (EventWriterIter = copy.mEventWriterCollection->begin();
       EventWriterIter != copy.mEventWriterCollection->end();
       EventWriterIter++) {
    mEventWriterCollection->push_back(*EventWriterIter);
  }
}

//_________________
MpdFemtoManager& MpdFemtoManager::operator=(const MpdFemtoManager& man) {
  // Assignment operator
  if (this != &man) {

    mEventReader = man.mEventReader;

    // Clean collections
    MpdFemtoAnalysisIterator analysisIter;
    if (mAnalysisCollection) {
      for (analysisIter = mAnalysisCollection->begin();
	   analysisIter != mAnalysisCollection->end(); analysisIter++) {
	delete *analysisIter;
	*analysisIter = nullptr;
      }
    } //if (mAnalysisCollection)
    delete mAnalysisCollection;

    MpdFemtoEventWriterIterator writerIter;
    if (mEventWriterCollection) {
      for (writerIter = mEventWriterCollection->begin();
	   writerIter != mEventWriterCollection->end(); writerIter++) {
	delete *writerIter;
	*writerIter = nullptr;
      }
    } //if(mEventWriterCollection)
    delete mEventWriterCollection;

    // Copy collections
    mAnalysisCollection = new MpdFemtoAnalysisCollection;
    for (analysisIter = man.mAnalysisCollection->begin();
	 analysisIter != man.mAnalysisCollection->end(); analysisIter++) {
      mAnalysisCollection->push_back(*analysisIter);
    }

    mEventWriterCollection = new MpdFemtoEventWriterCollection;
    for (writerIter = man.mEventWriterCollection->begin();
	 writerIter != man.mEventWriterCollection->end(); writerIter++) {
      mEventWriterCollection->push_back(*writerIter);
    }

  } //if ( this != &man )

  return *this;
}

//_________________
MpdFemtoManager::~MpdFemtoManager() {
  // Destructor
  delete mEventReader;

  // Delete each Analysis in the Collection
  // and then the Collection itself
  MpdFemtoAnalysisIterator AnalysisIter;
  for (AnalysisIter = mAnalysisCollection->begin();
       AnalysisIter != mAnalysisCollection->end();
       AnalysisIter++) {
    delete *AnalysisIter;
    *AnalysisIter = nullptr;
  }
  delete mAnalysisCollection;

  // Delete each EventWriter in the Collection,
  // and then the Collection itself
  MpdFemtoEventWriterIterator EventWriterIter;
  for (EventWriterIter = mEventWriterCollection->begin();
       EventWriterIter != mEventWriterCollection->end();
       EventWriterIter++) {
    delete *EventWriterIter;
    *EventWriterIter = nullptr;
  }
  delete mEventWriterCollection;
}

//_________________
int MpdFemtoManager::init() {
  MpdFemtoString readerMessage;
  readerMessage += "*** *** *** *** *** *** *** *** *** *** *** *** \n";
  // EventReader
  if (mEventReader) {
    if (mEventReader->init("r", readerMessage)) {
      std::cout << " MpdFemtoManager::init() - Reader initialization failed " << std::endl;
      return 1;
    }
    readerMessage += mEventReader->report();
  }

  // EventWriters
  MpdFemtoEventWriterIterator EventWriterIter;
  for (EventWriterIter = mEventWriterCollection->begin();
       EventWriterIter != mEventWriterCollection->end();
       EventWriterIter++) {
    // The message (MpdFemtoString) passed into Init will be at the file header.
    // for that reason take the readerReport, add my own report and pass as message
    MpdFemtoString writerMessage = readerMessage;
    writerMessage += "*** *** *** *** *** *** *** *** *** *** *** *** \n";
    writerMessage += (*EventWriterIter)->report();
    if (*EventWriterIter) {
      if ((*EventWriterIter)->init("w", writerMessage)) {
	// yes, the message from the reader is passed into the writer
	std::cout << " MpdFemtoManager::Init() - Writer initialization failed " << std::endl;
	return 1;
      } //if ( (*EventWriterIter)->Init("w",writerMessage) )
    } //if (*EventWriterIter)
  }
  return 0;
}

//_________________
void MpdFemtoManager::finish() {
  // EventReader
  if (mEventReader) {
    mEventReader->finish();
  }

  // EventWriters
  MpdFemtoEventWriterIterator EventWriterIter;
  MpdFemtoBaseEventWriter* currentEventWriter;
  for (EventWriterIter = mEventWriterCollection->begin();
       EventWriterIter != mEventWriterCollection->end();
       EventWriterIter++) {
    currentEventWriter = *EventWriterIter;
    currentEventWriter->finish();
  }

  // Analyses
  MpdFemtoAnalysisIterator AnalysisIter;
  MpdFemtoBaseAnalysis* currentAnalysis;
  for (AnalysisIter = mAnalysisCollection->begin();
       AnalysisIter != mAnalysisCollection->end();
       AnalysisIter++) {
    currentAnalysis = *AnalysisIter;
    currentAnalysis->finish();
  }
}

//_________________
MpdFemtoString MpdFemtoManager::report() {
  // Report construction
  string stemp;
  char ctemp[100];

  // EventReader
  stemp = mEventReader->report();

  // EventWriters
  sprintf(ctemp, "\nMpdFemtoManager Reporting %u EventWriters\n", (unsigned int) mEventWriterCollection->size());
  stemp += ctemp;
  MpdFemtoEventWriterIterator EventWriterIter;
  MpdFemtoBaseEventWriter* currentEventWriter;
  for (EventWriterIter = mEventWriterCollection->begin();
       EventWriterIter != mEventWriterCollection->end();
       EventWriterIter++) {
    std::cout << "MpdFemtoManager - asking for EventWriter Report" << std::endl;
    currentEventWriter = *EventWriterIter;
    stemp += currentEventWriter->report();
  }

  // Analyses
  sprintf(ctemp, "\nMpdFemtoManager Reporting %u Analyses\n", (unsigned int) mAnalysisCollection->size());
  stemp += ctemp;
  MpdFemtoAnalysisIterator AnalysisIter;
  MpdFemtoBaseAnalysis* currentAnalysis;
  for (AnalysisIter = mAnalysisCollection->begin();
       AnalysisIter != mAnalysisCollection->end();
       AnalysisIter++) {
    std::cout << "MpdFemtoManager - asking for Analysis Report" << std::endl;
    currentAnalysis = *AnalysisIter;
    stemp += currentAnalysis->report();
  }

  MpdFemtoString returnThis = stemp;
  return returnThis;
}

//_________________
MpdFemtoBaseAnalysis* MpdFemtoManager::analysis(int n) {
  // Return pointer to n-th analysis
  if (n < 0 || n > (int) mAnalysisCollection->size()) {
    return nullptr;
  }
  MpdFemtoAnalysisIterator iter = mAnalysisCollection->begin();
  for (int i = 0; i < n; i++) {
    iter++;
  }
  return *iter;
}

//_________________
MpdFemtoBaseEventWriter* MpdFemtoManager::eventWriter(int n) {
  // Return pointer to n-th analysis
  if (n < 0 || n > (int) mEventWriterCollection->size()) {
    return nullptr;
  }
  MpdFemtoEventWriterIterator iter = mEventWriterCollection->begin();
  for (int i = 0; i < n; i++) {
    iter++;
  }
  return *iter;
}

//_________________
int MpdFemtoManager::processEvent() {
  // Process a single event by reading it and passing it to each
  // analysis and event writer
  // NOTE - this ReturnHbtEvent makes a *new* MpdFemtoEvent - delete it when done!
  MpdFemtoEvent* currentHbtEvent = mEventReader->returnHbtEvent();

  // Uncomment next line for debug only!
  //  std::cout << "Event reader has returned control to manager" << std::endl;

  // If no HbtEvent is returned, then we abort processing.
  // the question is now: do we try again next time (i.e. there may be an HbtEvent next time)
  // or are we at EOF or something?  If Reader says Status=0, then that means try again later.
  // so, we just return the Reader's Status.

  if (!currentHbtEvent) {
    // Reader may return NULL pointer in case it doesn't receive file
    // or current event doesn't pass triggerId cut
    // So, in order not to flood output next line is commented

    // Uncomment for debug only!
    //  std::cout << "MpdFemtoManager::processEvent() - Reader::ReturnHbtEvent() has returned null pointer\n";

    return mEventReader->status();
  }

  // Loop over all the EventWriters
  MpdFemtoEventWriterIterator EventWriterIter;
  for (EventWriterIter = mEventWriterCollection->begin();
       EventWriterIter != mEventWriterCollection->end();
       EventWriterIter++) {

    // Uncomment next line for debug only
    //  std::cout << " *EventWriterIter " <<  *EventWriterIter << std::endl;

    (*EventWriterIter)->writeHbtEvent(currentHbtEvent);
  }

  // Loop over all the Analysis
  MpdFemtoAnalysisIterator AnalysisIter;
  for (AnalysisIter = mAnalysisCollection->begin();
       AnalysisIter != mAnalysisCollection->end();
       AnalysisIter++) {

    // Uncomment next line for debug only
    // std::cout << "MpdFemtoManager::processEvent() - Calling processEvent for analysis " << std::endl;
    (*AnalysisIter)->processEvent(currentHbtEvent);
  }

  if (currentHbtEvent) {
    // Uncomment next line for debug only
    // std::cout << "MpdFemtoManager::processEvent() - Removing event";
    delete currentHbtEvent;
    currentHbtEvent = nullptr;
    // Uncomment for debug only
    //std::cout << "\t[DONE]" << std::endl;
  }

  // Uncomment for debug only
  // std::cout << "MpdFemtoManager::processEvent() - return to caller ... " << std::endl;

  return 0; // 0 = "good return"
}
