//
// The MpdMiniArrays holds names, types and initial sizes of the mini arrays
//

// MiniDst headers
#include "MpdMiniArrays.h"

//    ARRAY NAMES
//_________________
const char *MpdMiniArrays::miniArrayNames[NAllMiniArrays] = {
   "Event", "Track", "BTofHit", "BTofPidTraits", "BECalCluster", "TrackCovMatrix", "FHCalHit", "McEvent", "McTrack"};

//   ARRAY TYPES
//_________________
const char *MpdMiniArrays::miniArrayTypes[NAllMiniArrays] = {
   "MpdMiniEvent",          "MpdMiniTrack",    "MpdMiniBTofHit", "MpdMiniBTofPidTraits", "MpdMiniBECalCluster",
   "MpdMiniTrackCovMatrix", "MpdMiniFHCalHit", "MpdMiniMcEvent", "MpdMiniMcTrack"};

//              ARRAY SIZES
// These are intial sizes. Automatically resized if too small.
// Choosing too large initial values gives a performance penalty when reading
// only selected miniDst branches
//_________________
int MpdMiniArrays::miniArraySizes[NAllMiniArrays] = {
   1,   // MpdMiniEvent
   800, // MpdMiniTrack
   100, // MpdMiniBTofHit
   100, // MpdMiniBTofPidTraits
   300, // MpdMiniBECalCluster
   800, // MpdMiniTrackCovMatrix
   616, // MpdMiniFHCalHit
   1,   // MpdMiniMcEvent
   1000 // MpdMiniMcTrack
};

//_________________
MpdMiniArrays::MpdMiniArrays()
{
   // Default constructor
   /* empty */
}

//_________________
MpdMiniArrays::~MpdMiniArrays()
{
   // Destructor
   /* empty */
}
