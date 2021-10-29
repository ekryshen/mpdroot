/**
 * \class MpdFemtoEventCutCollection
 * \brief Holds collection of event cuts
 *
 * The class holds collection (STL list) of event cuts
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoEventCutCollection_h
#define MpdFemtoEventCutCollection_h

// C++ headers
#include <list>

// Forward declaration
class MpdFemtoBaseEventCut;

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef std::list<MpdFemtoBaseEventCut*, std::allocator<MpdFemtoBaseEventCut*> > MpdFemtoEventCutCollection;
typedef std::list<MpdFemtoBaseEventCut*, std::allocator<MpdFemtoBaseEventCut*> >::iterator MpdFemtoEventCutIterator;
#else
typedef std::list<MpdFemtoBaseEventCut*> MpdFemtoEventCutCollection;
typedef std::list<MpdFemtoBaseEventCut*>::iterator MpdFemtoEventCutIterator;
#endif

#endif // #define MpdFemtoEventCutCollection_h
