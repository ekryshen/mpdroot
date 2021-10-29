/**
 * \class MpdFemtoTrackCollection
 * \brief Holds collection of tracks
 *
 * Holds a collectin (STL list) of tracks
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoTrackCollection_h
#define MpdFemtoTrackCollection_h

// C++ headers
#include <list>

// MpdFemtoMaker headers
#include "MpdFemtoTrack.h"

#if !defined(ST_NO_NAMESPACES)
using std::list;
#endif

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef list<MpdFemtoTrack*, allocator<MpdFemtoTrack*> > MpdFemtoTrackCollection;
typedef list<MpdFemtoTrack*, allocator<MpdFemtoTrack*> >::iterator MpdFemtoTrackIterator;
#else
typedef list<MpdFemtoTrack*> MpdFemtoTrackCollection;
typedef list<MpdFemtoTrack*>::iterator MpdFemtoTrackIterator;
#endif

#endif // #define MpdFemtoTrackCollection_h
