/**
 * \class MpdFemtoTrackCollection
 * \brief Holds collection of Xi(s)
 *
 * Holds a collectin (STL list) of Xi(s)
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoXiCollection_h
#define MpdFemtoXiCollection_h

// C++ headers
#include <list>

// MpdFemtoMaker headers
#include "MpdFemtoXi.h"

#if !defined(ST_NO_NAMESPACES)
using std::list;
#endif

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef list<MpdFemtoXi*, allocator<MpdFemtoXi*> > MpdFemtoXiCollection;
typedef list<MpdFemtoXi*, allocator<MpdFemtoXi*> >::iterator MpdFemtoXiIterator;
#else
typedef list<MpdFemtoXi*> MpdFemtoXiCollection;
typedef list<MpdFemtoXi*>::iterator MpdFemtoXiIterator;
#endif

#endif // #define MpdFemtoXiCollection_h
