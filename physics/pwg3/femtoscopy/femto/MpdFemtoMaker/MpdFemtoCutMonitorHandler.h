/**
 * \class MpdFemtoCutMonitorHandler
 * \brief A handler for cut monitors
 *
 * You add cut monitors to the collection which are stored in two separate
 * collections - one which stores characteristics of the entities (tracks,
 * particles, pairs, events) that pass the respective cuts and the other for
 * the ones that fail the cut.
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoCutMonitorHandler_h
#define MpdFemtoCutMonitorHandler_h

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseCutMonitor.h"

// Infrastructure
#include "MpdFemtoTypes.h"
#include "MpdFemtoEvent.h"
#include "MpdFemtoTrack.h"
#include "MpdFemtoV0.h"
#include "MpdFemtoKink.h"
#include "MpdFemtoXi.h"
#include "MpdFemtoPair.h"
#include "MpdFemtoParticleCollection.h"
#include "MpdFemtoCutMonitorCollection.h"

// ROOT headers
#include "Rtypes.h"
#include "TList.h"

//_________________

class MpdFemtoCutMonitorHandler {
public:
    /// Default constructor
    MpdFemtoCutMonitorHandler();
    /// Copy constructor
    MpdFemtoCutMonitorHandler(const MpdFemtoCutMonitorHandler& copy);
    /// Assignment operator
    MpdFemtoCutMonitorHandler& operator=(const MpdFemtoCutMonitorHandler& copy);
    /// Destructor
    virtual ~MpdFemtoCutMonitorHandler();

    /// Return monitor collection that passed cuts

    MpdFemtoCutMonitorCollection* passMonitorColl() {
        return mPassColl;
    }
    /// Return monitor collection that failed to pass cuts

    MpdFemtoCutMonitorCollection* failMonitorColl() {
        return mFailColl;
    }
    /// Return n-th monitor that passed cut
    MpdFemtoBaseCutMonitor* passMonitor(int n);
    /// Return n-th monitor that failed to pass cut
    MpdFemtoBaseCutMonitor* failMonitor(int n);

    /// Add cut monitor
    void addCutMonitor(MpdFemtoBaseCutMonitor* cutMoni1, MpdFemtoBaseCutMonitor* cutMoni2);
    /// Add cut monitor
    void addCutMonitor(MpdFemtoBaseCutMonitor* cutMoni);
    /// Add cut monitor that will be written in case the cut will be passed
    void addCutMonitorPass(MpdFemtoBaseCutMonitor* cutMoni);
    /// Add cut monitor that will be written in case the cut will not be passed
    void addCutMonitorFail(MpdFemtoBaseCutMonitor* cutMoni);

    /// Fill cut monitor for the event
    void fillCutMonitor(const MpdFemtoEvent* event, bool pass);
    /// Fill cut monitor for the track
    void fillCutMonitor(const MpdFemtoTrack* track, bool pass);
    /// Fill cut monitor for the v0
    void fillCutMonitor(const MpdFemtoV0* v0, bool pass);
    /// Fill cut monitor for the kink
    void fillCutMonitor(const MpdFemtoKink* kink, bool pass);
    /// Fill cut monitor for the xi
    void fillCutMonitor(const MpdFemtoXi* xi, bool pass);
    /// Fill cut monitor for the pair
    void fillCutMonitor(const MpdFemtoPair* pair, bool pass);
    /// Fill cut monitor for the pair
    void fillCutMonitor(const MpdFemtoParticleCollection* partColl);
    /// Fill cut monitor for the event and collection
    void fillCutMonitor(const MpdFemtoEvent* event, const MpdFemtoParticleCollection* partColl);
    /// Fill cut monitor for two collections
    void fillCutMonitor(const MpdFemtoParticleCollection* partColl1, const MpdFemtoParticleCollection* partColl2);

    /// Call finish
    void finish();

    /// Obtain list of objects to be written as an output
    virtual TList *getOutputList();

    /// Event begin
    virtual void eventBegin(const MpdFemtoEvent* event);
    /// Event end
    virtual void eventEnd(const MpdFemtoEvent* event);

private:

    /// Are the collections empty?
    bool mCollectionsEmpty;
    /// Collection of cut monitors for passed entities
    MpdFemtoCutMonitorCollection* mPassColl;
    /// Collection of cut monitors for failed entities
    MpdFemtoCutMonitorCollection* mFailColl;

    ClassDef(MpdFemtoCutMonitorHandler, 0)
};

#endif // #define MpdFemtoCutMonitorHandler_h
