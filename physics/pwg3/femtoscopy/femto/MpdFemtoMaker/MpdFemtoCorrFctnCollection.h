/**
 * \class MpdFemtoCorrFctnCollection
 * \brief Holds collection of correlation functions
 *
 *   The CorrFctnCollection contains pointers to all correlation functions
 *   that are associated with a particular analysis object.
 */

#ifndef MpdFemtoCorrFctnCollection_h
#define MpdFemtoCorrFctnCollection_h

// C++ headeres
#include <list>
#if !defined(ST_NO_NAMESPACES)
using std::list;
#endif

// Forward declaration
class MpdFemtoBaseCorrFctn;

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef list<MpdFemtoBaseCorrFctn*, allocator<MpdFemtoBaseCorrFctn*> > MpdFemtoCorrFctnCollection;
typedef list<MpdFemtoBaseCorrFctn*, allocator<MpdFemtoBaseCorrFctn*> >::iterator MpdFemtoCorrFctnIterator;
#else
typedef list<MpdFemtoBaseCorrFctn*> MpdFemtoCorrFctnCollection;
typedef list<MpdFemtoBaseCorrFctn*>::iterator MpdFemtoCorrFctnIterator;
#endif

#endif // #define MpdFemtoCorrFctnCollection_h
