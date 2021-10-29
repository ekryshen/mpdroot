/**
 * \class MpdFemtoKinkCollection
 * \brief Holds collection of kinks
 *
 * Holds a collectin (STL list) of kinks
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoKinkCollection_h
#define MpdFemtoKinkCollection_h

// C++ headers
#include <list>

// MpdFemtoMaker headers
#include "MpdFemtoKink.h"

#if !defined(ST_NO_NAMESPACES)
using std::list;
#endif

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef list<MpdFemtoKink*, allocator<MpdFemtoKink*> > MpdFemtoKinkCollection;
typedef list<MpdFemtoKink*, allocator<MpdFemtoKink*> >::iterator MpdFemtoKinkIterator;
#else
typedef list<MpdFemtoKink*> MpdFemtoKinkCollection;
typedef list<MpdFemtoKink*>::iterator MpdFemtoKinkIterator;
#endif

#endif // #define MpdFemtoKinkCollection_h
