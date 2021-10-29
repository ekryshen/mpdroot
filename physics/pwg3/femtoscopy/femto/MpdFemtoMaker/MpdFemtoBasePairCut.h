/**
 * \class MpdFemtoBasePairCut
 * \brief The pure virtual base class for the pair cut
 *
 * All pair cuts must inherit from this one
 *
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoBasePairCut_h
#define MpdFemtoBasePairCut_h

// C++ headers
#include <string>

// Forward declaration
class MpdFemtoBaseAnalysis;

// MpdFemtoMaker headers
#include "MpdFemtoString.h"
#include "MpdFemtoEvent.h"
#include "MpdFemtoPair.h"
#include "MpdFemtoCutMonitorHandler.h"

// ROOT headers
#include "TList.h"
#include "TObjString.h"

//_________________

class MpdFemtoBasePairCut : public MpdFemtoCutMonitorHandler {
public:
    /// Default constructor
    MpdFemtoBasePairCut();
    /// Copy constructor
    MpdFemtoBasePairCut(const MpdFemtoBasePairCut& c);
    /// Assignment operator
    MpdFemtoBasePairCut& operator=(const MpdFemtoBasePairCut& c);
    /// Default destructor

    virtual ~MpdFemtoBasePairCut() {
        /* no-op */
    }

    /// Returns true in pair passed the cut and false if not
    virtual bool pass(const MpdFemtoPair* pair) = 0;

    /// User-written method to return string describing cuts
    virtual MpdFemtoString report() = 0;
    /// Returns a TList with settings
    virtual TList *listSettings() = 0;

    /// Declare event start

    virtual void eventBegin(const MpdFemtoEvent*) {
        /* no-op */
    }
    /// Declare event end

    virtual void eventEnd(const MpdFemtoEvent*) {
        /* no-op */
    }
    /// Clone pair cut

    virtual MpdFemtoBasePairCut* clone() {
        return nullptr;
    }

    /// The following allows "back-pointing" from the CorrFctn
    /// to the "parent" Analysis
    friend class MpdFemtoBaseAnalysis;
    /// Return pointer to the analysis

    MpdFemtoBaseAnalysis* hbtAnalysis() {
        return mBaseAnalysis;
    }
    /// Set analysis

    void setAnalysis(MpdFemtoBaseAnalysis* ana) {
        mBaseAnalysis = ana;
    }

protected:
    /// Pointer to the base analysis
    MpdFemtoBaseAnalysis* mBaseAnalysis; //!<!

    ClassDef(MpdFemtoBasePairCut, 0)
};

//_________________

inline MpdFemtoBasePairCut::MpdFemtoBasePairCut() : MpdFemtoCutMonitorHandler(), mBaseAnalysis() {
    /* empty */
}

//_________________

inline MpdFemtoBasePairCut::MpdFemtoBasePairCut(const MpdFemtoBasePairCut& /* c */) : MpdFemtoCutMonitorHandler(), mBaseAnalysis(nullptr) {
    /* empty */
}

//_________________

inline MpdFemtoBasePairCut& MpdFemtoBasePairCut::operator=(const MpdFemtoBasePairCut& c) {
    if (this != &c) {
        mBaseAnalysis = c.mBaseAnalysis;
    }
    return *this;
}

#endif // #define MpdFemtoBasePairCut_h
