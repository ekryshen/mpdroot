/**
 * \class MpdMiniArray
 * \brief Stores mini arrays
 *
 * MpdMiniArrays is a pure C++ class that holds names of the mini arrays,
 * names of TBranches and TClones arrays.
 *
 * The currently implemented objects:
 * Event          - Event
 * Track          - Track (reconstructed for both real and MC)
 * BTofHit        - Hit in Barrel Time-Of-Flight
 * BTofPidTraits  - PID information from TOF for the TOF-matched track
 * BECalHit       - Hit in Barrel Electromagnectic Calorimeter (for each tower)
 * BECalPidTraits - PID information from the tower for the ECal-matched track
 * TrackCovMatrix - Track covariance matrix
 * FHCalHit       - Hit in forward hadronic calorimeter (for each tower)
 * McEvent        - Monte Carlo event information
 * McTrack        - Monte Carlo track information
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \email nigmatkulov@gmail.com ; ganigmatkulov@mephi.ru
 * \date May 01, 2020
 */

#ifndef MpdMiniArrays_h
#define MpdMiniArrays_h

//_________________
class MpdMiniArrays {

public:
   /// Default constructor
   MpdMiniArrays();

   /// Destructor
   virtual ~MpdMiniArrays();

   /// Should be changed to constexpr once ROOT 6 is available at STAR
   enum { NAllMiniArrays = 9 };

   /// Names of the TBranches in the TTree/File
   static const char *miniArrayNames[NAllMiniArrays];

   /// Names of the classes, the TClonesArrays are arrays of this type
   static const char *miniArrayTypes[NAllMiniArrays];

   /// Maximum sizes of the TClonesArrays
   static int miniArraySizes[NAllMiniArrays];

   /// Array names
   enum TypeIndex {
      Event = 0,
      Track,
      BTofHit,
      BTofPidTraits,
      BECalCluster,
      TrackCovMatrix,
      FHCalHit,
      McEvent,
      McTrack
   };
};

#endif // #define MpdMiniArrays_h
