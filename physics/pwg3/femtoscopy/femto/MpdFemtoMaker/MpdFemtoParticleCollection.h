/**
 * \class MpdFemtoParticleCollection
 * \brief Holds collection of particles
 *
 * Holds collection (STL list) of particles
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoParticleCollection_h
#define MpdFemtoParticleCollection_h

// C++ headers
#include <list>

// MpdFemtoMaker headers
#include "MpdFemtoParticle.h"

#if !defined(ST_NO_NAMESPACES)
using std::list;
#endif

#ifdef ST_NO_TEMPLATE_DEF_ARGS
typedef list<MpdFemtoParticle*, allocator<MpdFemtoParticle*> > MpdFemtoParticleCollection;
typedef list<MpdFemtoParticle*, allocator<MpdFemtoParticle*> >::iterator MpdFemtoParticleIterator;
typedef list<MpdFemtoParticle*, allocator<MpdFemtoParticle*> >::const_iterator MpdFemtoParticleConstIterator;
#else
typedef list<MpdFemtoParticle*> MpdFemtoParticleCollection;
typedef list<MpdFemtoParticle*>::iterator MpdFemtoParticleIterator;
typedef list<MpdFemtoParticle*>::const_iterator MpdFemtoParticleConstIterator;
#endif

#endif // #define MpdFemtoParticleCollection_h
