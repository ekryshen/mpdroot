/**
 * \class MpdFemtoV0Collection
 * \brief Holds collection of V0s
 *
 * Holds a collectin (STL list) of V0s
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoV0Collection_h
#define MpdFemtoV0Collection_h

// C++ headers
#include <list>

// MpdFemtoMaker headers
#include "MpdFemtoV0.h"

#if !defined(ST_NO_NAMESPACES)
using std::list;
#endif

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef list<MpdFemtoV0*, allocator<MpdFemtoV0*> > MpdFemtoV0Collection;
typedef list<MpdFemtoV0*, allocator<MpdFemtoV0*> >::iterator MpdFemtoV0Iterator;
#else
typedef list<MpdFemtoV0*> MpdFemtoV0Collection;
typedef list<MpdFemtoV0*>::iterator MpdFemtoV0Iterator;
#endif

#endif // #define MpdFemtoV0Collection_h
