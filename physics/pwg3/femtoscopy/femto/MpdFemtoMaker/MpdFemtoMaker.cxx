//
// The Maker that runs MpdFemtoMaker maker subclasses.
//

// C++ headers
#include <iostream>
#include <stdio.h>
#include <string>

// ROOT headers
#include "Riostream.h"

// MpdFemtoMaker headers
#include "MpdFemtoMaker.h"

ClassImp(MpdFemtoMaker)

//_________________
MpdFemtoMaker::MpdFemtoMaker(const char* name __attribute__ ((unused)),
			     const char* title __attribute__ ((unused))) :
#ifdef sl73_gcc485
// STAR compiled library
StMaker(name, title),
#endif
mDebug(0) {
  // Constructor
  mHbtManager = new MpdFemtoManager;
  std::cout << string::npos << std::endl;
}

//_________________
MpdFemtoMaker::~MpdFemtoMaker() {
  // Destructor
  std::cout << "Inside ReaderMaker Destructor" << std::endl;
#ifdef sl73_gcc485
  SafeDelete(mHbtManager); //
#else
  delete mHbtManager;
#endif
}

//_________________
void MpdFemtoMaker::clear(const char*) {
  /* no-op - do not delete manager! */
#ifdef sl73_gcc485
  StMaker::Clear();
#endif
}

//_________________
Int_t MpdFemtoMaker::init() {
  if (mHbtManager->init()) {
    std::cout << "MpdFemtoMaker::Init() - manager init failed " << std::endl;
    return (2);
  }

  std::cout << "MpdFemtoMaker::Init() - requesting a report " << std::endl;

  MpdFemtoString tempString = mHbtManager->report();
  std::cout << "Got the report, now let me try to put it to screen" << std::endl;
  std::cout << tempString.c_str() << std::endl; //!

#ifdef sl73_gcc485
  return StMaker::Init();
#else
  return 0;
#endif
}

//_________________
Int_t MpdFemtoMaker::finish() {
  std::cout << mHbtManager->report().c_str() << std::endl;
  mHbtManager->finish();
#ifdef sl73_gcc485
  return StMaker::Finish();
#else
  return 0;
#endif
}

//_________________
Int_t MpdFemtoMaker::make() {
  if (mDebug > 1) {
    std::cout << "\nMpdFemtoMaker::Make -- processing event" << std::endl;
  }
#ifdef sl73_gcc485
  if (mHbtManager->processEvent()) {
    return kStEOF; // non-zero value returned --> end of file action
  } else {
    return kStOK;
  }
#else
  return mHbtManager->processEvent();
#endif
}
