/**
 * \class MpdFemtoMaker
 * \brief The Maker that runs MpdFemtoMaker maker subclasses.
 *
 * One uses the maker to instantiate and add analysis into it.
 * To work in the STAR environment, sl73_gcc485 global variable
 * must be defined at the compilation and a the run time.
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoMaker_h
#define MpdFemtoMaker_h

#ifdef sl73_gcc485
// STAR headers
#include "StMaker.h"
#endif

// MpdFemtoMaker headers
#include "MpdFemtoManager.h"

//_________________
class MpdFemtoMaker
#ifdef sl73_gcc485
: public StMaker
#endif
{
 public:
  /// Constructor
  MpdFemtoMaker(const char* name = "MpdFemto", const char* title = "MpdFemtoTitle");
  /// Destructor
  virtual ~MpdFemtoMaker();

  // Standard StChain (STAR) options that can be used
  // in a standalone mode.

  /// Clear method
  virtual void clear(const char* opt = "");
  /// Clear method
  virtual void Clear(const char* opt) {
    clear(opt);
  }
  /// Init method
  virtual Int_t init();
  /// Init method
  virtual Int_t Init() {
    return init();
  }
  /// Make method
  virtual Int_t make();
  /// Make method
  virtual Int_t Make() {
    return make();
  }
  /// Finish method
  virtual Int_t finish();
  /// Finish method

  virtual Int_t Finish() {
    return finish();
  }

  /// Return a pointer to MpdFemtoManager
  MpdFemtoManager* hbtManager() {
    return mHbtManager;
  }
  /// Return a pointer to MpdFemtoManager
  MpdFemtoManager* HbtManager() {
    return hbtManager();
  }
  /// Return debug value
  int debug() const {
    return mDebug;
  }
  /// Return debug value
  int Debug() const {
    return debug();
  }
  /// Set debug value
  void setDebug(int debug) {
    mDebug = debug;
  }
  /// Set debug value
  void SetDebug(int debug) {
    setDebug(debug);
  }

#ifdef sl73_gcc485
  /// Pointer to a chain
  StMaker* currentChain;
#endif

 private:

  /// Pointer to MpdFemtoManager
  MpdFemtoManager* mHbtManager; //! tells cint to skip it

  /// Debug value
  int mDebug;

  ClassDef(MpdFemtoMaker, 0)
};

#endif // #define MpdFemtoMaker_h
