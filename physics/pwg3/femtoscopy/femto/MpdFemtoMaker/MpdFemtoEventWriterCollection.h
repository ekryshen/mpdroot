/**
 * \class MpdFemtoEventWriterCollection
 * \brief Holds collection of Writers
 *
 * The class keeps a collection (STL list) of event writers
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoEventWriterCollection_h
#define MpdFemtoEventWriterCollection_h

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseEventWriter.h"

// C++ headers
#include <list>

#if !defined(ST_NO_NAMESPACES)
using std::list;
#endif

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef list<MpdFemtoBaseEventWriter*, allocator<MpdFemtoBaseEventWriter*> > MpdFemtoEventWriterCollection;
typedef list<MpdFemtoBaseEventWriter*, allocator<MpdFemtoBaseEventWriter*> >::iterator MpdFemtoEventWriterIterator;
#else
typedef list<MpdFemtoBaseEventWriter*> MpdFemtoEventWriterCollection;
typedef list<MpdFemtoBaseEventWriter*>::iterator MpdFemtoEventWriterIterator;
#endif

#endif // #define MpdFemtoEventWriterCollection_h
