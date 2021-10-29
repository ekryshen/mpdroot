/**
 * \class MpdFemtoBaseTrackCut
 * \brief Base class for track cuts
 *
 * The pure virtual base class for track cuts. All track cuts must
 * inherit from this one.
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoBaseTrackCut_h
#define MpdFemtoBaseTrackCut_h

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseParticleCut.h"
// Infrastructure
#include "MpdFemtoTypes.h"
#include "MpdFemtoTrack.h"

//_________________

class MpdFemtoBaseTrackCut : public MpdFemtoBaseParticleCut {
public:
    /// Default constructor
    MpdFemtoBaseTrackCut();
    /// Copy constructor
    MpdFemtoBaseTrackCut(const MpdFemtoBaseTrackCut& copy);
    /// Assignment operator
    MpdFemtoBaseTrackCut& operator=(const MpdFemtoBaseTrackCut& copy);
    /// Default destructor

    virtual ~MpdFemtoBaseTrackCut() {
        /* no-op */
    }

    /// Returns true if passed the track cut and false if not
    virtual bool pass(const MpdFemtoTrack* track) = 0;
    /// Return track type

    virtual MpdFemtoParticleType type() {
        return hbtTrack;
    }
    /// Clone track cut

    virtual MpdFemtoBaseTrackCut* clone() {
        return nullptr;
    }

    ClassDef(MpdFemtoBaseTrackCut, 0)
};

//_________________

inline MpdFemtoBaseTrackCut::MpdFemtoBaseTrackCut() : MpdFemtoBaseParticleCut() {
    /* empty */
}

//_________________

inline MpdFemtoBaseTrackCut::MpdFemtoBaseTrackCut(const MpdFemtoBaseTrackCut& c) : MpdFemtoBaseParticleCut(c) {
    /* empty */
}

//_________________

inline MpdFemtoBaseTrackCut& MpdFemtoBaseTrackCut::operator=(const MpdFemtoBaseTrackCut& c) {
    if (this != &c) {
        MpdFemtoBaseParticleCut::operator=(c);
    }
    return *this;
}

#endif // #define MpdFemtoTrackCut_h
