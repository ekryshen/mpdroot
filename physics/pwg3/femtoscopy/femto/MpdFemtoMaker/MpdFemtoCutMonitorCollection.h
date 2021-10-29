/**
 * \class MpdFemtoCutMonitorCollection
 * \brief Keeps collection of cut monitors
 *
 * The class keeps a collection of cut monitors stored as
 * STL vector
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoCutMonitorCollection_h
#define MpdFemtoCutMonitorCollection_h

// C++ headers
#include <vector>
#include <iterator>
//#include <memory> // for shared and unique pointers (should we use them?)

#if !defined(ST_NO_NAMESPACES)
using std::vector;
#endif

// MpdFemtoMaker headers
class MpdFemtoBaseCutMonitor;

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef vector<MpdFemtoBaseCutMonitor*, allocator<MpdFemtoBaseCutMonitor*> > MpdFemtoCutMonitorCollection;
typedef vector<MpdFemtoBaseCutMonitor*, allocator<MpdFemtoBaseCutMonitor*> >::iterator MpdFemtoCutMonitorIterator;
#else
typedef vector<MpdFemtoBaseCutMonitor*> MpdFemtoCutMonitorCollection;
typedef vector<MpdFemtoBaseCutMonitor*>::iterator MpdFemtoCutMonitorIterator;
#endif

#endif // #define MpdFemtoCutMonitorCollection_h
