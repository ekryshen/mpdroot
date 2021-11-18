set(INCLUDE_DIRECTORIES
${BASE_INCLUDE_DIRECTORIES}
${CMAKE_SOURCE_DIR}/mpddst/MpdMiniEvent
)

Set(SYSTEM_INCLUDE_DIRECTORIES
 ${ROOT_INCLUDE_DIR}
)

include_directories(${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
 ${ROOT_LIBRARY_DIR}
 ${FAIRROOT_LIBRARY_DIR}
)
 
link_directories(${LINK_DIRECTORIES})

# List of source files
set(SRCS
MpdMiniPhysicalHelix.cxx
MpdMiniEvent.cxx
MpdMiniDst.cxx
MpdMiniBECalCluster.cxx
MpdMiniArrays.cxx
MpdMiniTrack.cxx
MpdMiniTrackCovMatrix.cxx
MpdMiniMcTrack.cxx
MpdMiniMcEvent.cxx
MpdMiniHelix.cxx
MpdMiniDstReader.cxx
MpdMiniBTofPidTraits.cxx
MpdMiniBTofHit.cxx
MpdMiniFHCalHit.cxx
)

Set(HEADERS)
Set(LINKDEF MpdMiniEventLinkDef.h)
Set(LIBRARY_NAME MpdMiniEvent)
Set(DEPENDENCIES Base)

GENERATE_LIBRARY()
