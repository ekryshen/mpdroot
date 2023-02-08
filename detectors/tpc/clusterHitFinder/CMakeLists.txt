# Create a library called "libtpcClusterHitFinder"

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_CURRENT_LIST_DIR}
 ${CMAKE_CURRENT_LIST_DIR}/fast
 ${CMAKE_CURRENT_LIST_DIR}/mlem
 ${CMAKE_SOURCE_DIR}/detectors/tpc
 ${CMAKE_SOURCE_DIR}/detectors/tpc/digitizer
 ${CMAKE_SOURCE_DIR}/detectors/tpc/digitizer/AZlt
 ${CMAKE_SOURCE_DIR}/detectors/tpc/fairTpc
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
  fast/TpcClusterHitFinderFast.cxx
  hitProducer/MpdTpcHitProducer.cxx
  mlem/TpcClusterHitFinderMlem.cxx
  mlem/MpdTpc2dCluster.cxx
  mlem/MpdTpcHit.cxx
)

Set(HEADERS
  fast/TpcClusterHitFinderFast.h
  fast/TpcClustering.h
  hitProducer/MpdTpcHitProducer.h
  mlem/TpcClusterHitFinderMlem.h
  mlem/MpdTpc2dCluster.h
  mlem/MpdTpcHit.h
)

Set(LINKDEF tpcClusterHitFinderLinkDef.h)
Set(LIBRARY_NAME tpcClusterHitFinder)
Set(DEPENDENCIES Kalman MpdField MpdBase tpcGeometry)

GENERATE_LIBRARY()
