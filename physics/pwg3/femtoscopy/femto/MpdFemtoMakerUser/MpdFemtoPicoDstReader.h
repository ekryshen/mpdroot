/**
 * \class MpdFemtoPicoDstReader
 * \brief Reader for the StPicoDst format
 *
 * The reader class for StPicoDst format. It reads McDst and
 * converts data to the internal MpdFemtoMaker structure (MpdFemtoEvent,
 * MpdFemtoTrack, etc).
 */

#ifndef MpdFemtoPicoDstReader_h
#define MpdFemtoPicoDstReader_h

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseEventReader.h"
// Infrastructure
#include "MpdFemtoEvent.h"
#include "MpdFemtoTrack.h"
#include "MpdFemtoV0.h"
#include "MpdFemtoKink.h"
#include "MpdFemtoXi.h"
#include "MpdFemtoEnumeration.h"
#include "MpdFemtoString.h"

// StPicoDst headers
#include "StPicoDst.h"
#include "StPicoDstReader.h"
#include "StPicoEvent.h"
#include "StPicoTrack.h"
#include "StPicoBTofPidTraits.h"

// ROOT haders
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"

#ifndef _VANILLA_ROOT_
// Forward declarations
class StRefMultCorr;
#endif

//_________________

class MpdFemtoPicoDstReader : public MpdFemtoBaseEventReader {
public:
    /// Default constructor
    MpdFemtoPicoDstReader();
    /// Constructor
    MpdFemtoPicoDstReader(StPicoDstReader *picoDstReader, int debug = 1);
    /// Copy constructor
    MpdFemtoPicoDstReader(const MpdFemtoPicoDstReader& copy);
    /// Assignment operator
    MpdFemtoPicoDstReader& operator=(const MpdFemtoPicoDstReader& copy);
    /// Destructor
    virtual ~MpdFemtoPicoDstReader();

    /// Return MpdFemtoEvent (back to the analysis)
    virtual MpdFemtoEvent *returnHbtEvent();
    /// Make report
    virtual MpdFemtoString report();

    /// Set picoDst reader

    void setPicoDstReader(StPicoDstReader *picoDstReader) {
        mPicoDstReader = picoDstReader;
    }
#ifndef _VANILLA_ROOT_
    /// Set StRefMultCorr

    void setRefMultCorrUtil(StRefMultCorr *refMultCorr) {
        mRefMultCorrUtil = refMultCorr;
    }
#endif
    /// Set use StRefMultCorr package

    void setUseRefMultCorr(const bool& useRefMultCorr) {
        mUseRefMultCorr = useRefMultCorr;
    }
    /// Set use ZDC correction for StRefMultCorr package

    void setUseZdcCorrection(const bool& useZdcCorrection) {
        mUseZdcCorrection = useZdcCorrection;
    }

private:

    /// Sphericity calculation for tracks within |eta|<0.5 (or |eta|<1.),
    /// pT>0.15, DCA<3 and nHits>10
    void sphericityCalculation();

    /// Pointer to the StPicoDstReader
    StPicoDstReader *mPicoDstReader;
    /// Pointer to MpdFemtoEvent
    MpdFemtoEvent *mHbtEvent;

    /// Use StRefMultCorr
    bool mUseRefMultCorr;
    /// Use ZDC correction for StRefMultCorr
    bool mUseZdcCorrection;

    /// Sphericity esitmated by tracks with  |eta|<0.5, pT>0.15, DCA<3 and nHits>10
    float mSphericity;
    /// Sphericity esitmated by tracks with  |eta|<1, pT>0.15, DCA<3 and nHits>10
    float mSphericity2;

    /// Number of events passed
    long int mEventsPassed;

#ifndef _VANILLA_ROOT_
    /// Pointer to the StRefMultCorr instance
    StRefMultCorr *mRefMultCorrUtil;
#endif

    ClassDef(MpdFemtoPicoDstReader, 1)

};

#endif // #define MpdFemtoPicoDstReader_h
