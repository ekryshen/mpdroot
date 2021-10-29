/**
 * \class MpdFemtoPairCutCollection
 * \brief Holds collection of pair cuts
 *
 * Holds collection (STL list) of pair cuts
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoPairCutCollection_h
#define MpdFemtoPairCutCollection_h

// C++ headers
#include <list>

// Forward declarations
class MpdFemtoPairCut;

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef std::list<MpdFemtoPairCut*, std::allocator<MpdFemtoPairCut*> > MpdFemtoPairCutCollection;
typedef std::list<MpdFemtoPairCut*, std::allocator<MpdFemtoPairCut*> >::iterator MpdFemtoPairCutIterator;
#else
typedef std::list<MpdFemtoPairCut*> MpdFemtoPairCutCollection;
typedef std::list<MpdFemtoPairCut*>::iterator MpdFemtoPairCutIterator;
#endif

#endif // #define MpdFemtoPairCutCollection_h
