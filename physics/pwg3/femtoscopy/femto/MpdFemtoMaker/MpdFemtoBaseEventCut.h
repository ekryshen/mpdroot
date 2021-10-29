/**
 * \class MpdFemtoBaseEventCut
 * \brief The pure virtual base class for the event cut
 *
 * All event cuts must inherit from this one and implement the pass() and
 * report() methods. The ::clone() function simply returns nullptr, so if
 * users want their cuts to behave as expected, they should also write
 * their own.
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoBaseEventCut_h
#define MpdFemtoBaseEventCut_h

// Forward declarations
class MpdFemtoEvent;
class MpdFemtoBaseAnalysis;

// MpdFemtoMaker headers
#include "MpdFemtoCutMonitorHandler.h"
#include "MpdFemtoString.h"

// ROOT headers
#include "TList.h"
#include "TString.h"

//_________________
class MpdFemtoBaseEventCut : public MpdFemtoCutMonitorHandler {
 public:
  /// Default constructor
  MpdFemtoBaseEventCut();
  /// Copy constructor
  MpdFemtoBaseEventCut(const MpdFemtoBaseEventCut& copy);
  /// Assignment operator
  MpdFemtoBaseEventCut& operator=(const MpdFemtoBaseEventCut& copy);
  /// Default destructor
  virtual ~MpdFemtoBaseEventCut() { /* no-op */ }

  /// True if event has passed the cut and false if not
  virtual bool pass(const MpdFemtoEvent* event) = 0;

  /// Return new settings list.
  ///
  /// This method creates a new list of TObjStrings describing cut parameters.
  /// The default implementation automatically calls the AppendSettings method
  /// to fill the list, so users only need to overload that method.
  virtual TList* listSettings() const;

  /// Appends cut settings to a TList
  ///
  /// This method should be overloaded by the user to add any relevent settings
  /// of the cut to the list
  ///
  /// No settings are added by this class. Simply returns the incoming TList.
  ///
  /// \param A list to append settings to.
  /// \param prefix An optional prefix to prepend to the beginning of each setting
  /// \return The same pointer as the parameter
  virtual TList* appendSettings(TList*, const TString& prefix = "") const;

  /// User-written method to return string describing cuts
  virtual MpdFemtoString report() = 0;

  /// Default clone
  virtual MpdFemtoBaseEventCut* clone() const = 0;

  /// The following allows "back-pointing" from the CorrFctn
  ///to the "parent" Analysis
  friend class MpdFemtoBaseAnalysis;
  
  /// Return a pointer to the analysis
  MpdFemtoBaseAnalysis* hbtAnalysis() {
    return mBaseAnalysis;
  }
  
  /// Set analysis
  void setAnalysis(MpdFemtoBaseAnalysis* a) {
    mBaseAnalysis = a;
  }

 protected:
  /// Pointer to the base analysis
  MpdFemtoBaseAnalysis* mBaseAnalysis; //!<!

  ClassDef(MpdFemtoBaseEventCut, 0)
};

//_________________
inline MpdFemtoBaseEventCut::MpdFemtoBaseEventCut() : MpdFemtoCutMonitorHandler(), mBaseAnalysis(nullptr) {
  /* empty */
}

//_________________
inline MpdFemtoBaseEventCut::MpdFemtoBaseEventCut(const MpdFemtoBaseEventCut& c) :
			    MpdFemtoCutMonitorHandler(c), mBaseAnalysis(c.mBaseAnalysis) {
  /* empty */
}

//_________________
inline MpdFemtoBaseEventCut& MpdFemtoBaseEventCut::operator=(const MpdFemtoBaseEventCut& c) {
  if (this != &c) {
    mBaseAnalysis = c.mBaseAnalysis;
  }
  return *this;
}

//_________________
inline TList* MpdFemtoBaseEventCut::listSettings() const {
  return appendSettings(new TList());
}

//_________________
inline TList* MpdFemtoBaseEventCut::appendSettings(TList *listOfSettings, const TString& /* prefix */) const {
  return listOfSettings;
}

#endif // #define MpdFemtoBaseEventCut_h
