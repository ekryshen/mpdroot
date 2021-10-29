/**
 * \class MpdFemtoEnumeration
 * \brief Enumerations for MpdFemtoMaker
 *
 * The file keeps enumarations for particle types and I/O mode
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoEnumeration_h
#define MpdFemtoEnumeration_h

/// Particle types

enum MpdFemtoParticleType {
    hbtUndefined, hbtTrack, hbtV0, hbtKink, hbtXi
};
/// I/O mode

enum MpdFemtoIOMode {
    hbtRead, hbtWrite
};

#endif //#define MpdFemtoEnumeration_h
