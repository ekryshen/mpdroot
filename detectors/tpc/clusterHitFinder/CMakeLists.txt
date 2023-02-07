# Create a library called "libtpcClusterHitFinder"

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_CURRENT_LIST_DIR}
 ${CMAKE_CURRENT_LIST_DIR}/fast
 ${CMAKE_CURRENT_LIST_DIR}/mlem
 ${CMAKE_SOURCE_DIR}/detectors/tpc
 ${CMAKE_SOURCE_DIR}/detectors/tpc/digitizer
 ${CMAKE_SOURCE_DIR}/detectors/tpc/digitizer/AZlt
 ${CMAKE_SOURCE_DIR}/detectors/tpc/geometry
 ${CMAKE_SOURCE_DIR}/reconstruction/tracking/kalman
 ${CMAKE_SOURCE_DIR}/simulation/generators/mpdGen
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
  AbstractTpcClusterHitFinder.cxx
  fast/TpcClusterHitFinderFast.cxx
  mlem/TpcClusterHitFinderMlem.cxx
)

Set(HEADERS
  AbstractTpc2dCluster.h
  AbstractTpcClusterHitFinder.h
  AbstractTpcHit.h
  fast/TpcClusterHitFinderFast.h
  fast/TpcClustering.h
  mlem/TpcClusterHitFinderMlem.h
)

Set(LINKDEF tpcClusterHitFinderLinkDef.h)
Set(LIBRARY_NAME tpcClusterHitFinder)
Set(DEPENDENCIES MpdField MpdBase)

GENERATE_LIBRARY()
