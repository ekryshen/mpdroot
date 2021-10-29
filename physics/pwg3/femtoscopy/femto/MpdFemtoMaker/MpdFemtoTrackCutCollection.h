/**
 * \class MpdFemtoTrackCutCollection
 * \brief Holds collection of track cuts
 *
 * Holds collection (STL list) of track cuts
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoTrackCutCollection_h
#define MpdFemtoTrackCutCollection_h

// C++ headers
#include <list>

// MpdFemtoMaker headers
class MpdFemtoBaseTrackCut;

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef std::list<MpdFemtoBaseTrackCut*, std::allocator<MpdFemtoBaseTrackCut*> > MpdFemtoTrackCutCollection;
typedef std::list<MpdFemtoBaseTrackCut*, std::allocator<MpdFemtoBaseTrackCut*> >::iterator MpdFemtoTrackCutIterator;
#else
typedef std::list<MpdFemtoBaseTrackCut*> MpdFemtoTrackCutCollection;
typedef std::list<MpdFemtoBaseTrackCut*>::iterator MpdFemtoTrackCutIterator;
#endif

#endif // #define MpdFemtoTrackCutCollection_h
