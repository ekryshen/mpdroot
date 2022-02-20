set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_SOURCE_DIR}/mpddst/MpdMiniEvent
)

set(SYSTEM_INCLUDE_DIRECTORIES
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
  MpdMiniArrays.cxx
  MpdMiniBECalCluster.cxx
  MpdMiniBTofHit.cxx
  MpdMiniBTofPidTraits.cxx
  MpdMiniDst.cxx
  MpdMiniDstReader.cxx
  MpdMiniEvent.cxx
  MpdMiniFHCalHit.cxx
  MpdMiniHelix.cxx
  MpdMiniMcEvent.cxx
  MpdMiniMcTrack.cxx
  MpdMiniTrack.cxx
  MpdMiniTrackCovMatrix.cxx
  MpdMiniPhysicalHelix.cxx
)

set(HEADERS
  MpdMiniArrays.h
  MpdMiniBECalCluster.h
  MpdMiniBTofHit.h
  MpdMiniBTofPidTraits.h
  MpdMiniDst.h
  MpdMiniDstReader.h
  MpdMiniEvent.h
  MpdMiniFHCalHit.h
  MpdMiniHelix.h
  MpdMiniMcEvent.h
  MpdMiniMcTrack.h
  MpdMiniMessMgr.h
  MpdMiniPhysicalHelix.h
  MpdMiniTrack.h
  MpdMiniTrackCovMatrix.h
  PhysicalConstants.h
  SystemOfUnits.h
)

set(LINKDEF MpdMiniEventLinkDef.h)
set(LIBRARY_NAME MpdMiniEvent)
set(DEPENDENCIES Base)

GENERATE_LIBRARY()
