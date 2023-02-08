# Create a library called "libtpc" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_CURRENT_LIST_DIR}/clusterHitFinder
  ${CMAKE_CURRENT_LIST_DIR}/clusterHitFinder/fast
  ${CMAKE_CURRENT_LIST_DIR}/clusterHitFinder/mlem
  ${CMAKE_CURRENT_LIST_DIR}/digitizer
  ${CMAKE_CURRENT_LIST_DIR}/fairTpc
  ${CMAKE_CURRENT_LIST_DIR}/geometry
  ${CMAKE_SOURCE_DIR}/core/mpdField
  ${CMAKE_SOURCE_DIR}/detectors/tof
  ${CMAKE_SOURCE_DIR}/reconstruction/tracking/kalman
  ${CMAKE_SOURCE_DIR}/simulation/generators/mpdGen
  ${CMAKE_SOURCE_DIR}/simulation/mcStack
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
  MpdParticleIdentification.cxx
  MpdTpc2dCluster.cxx
  MpdTpcClusterFinderMlem.cxx
  MpdTpcDigit.cxx
  MpdTpcDigitizerAZlt.cxx
  MpdTpcDigitizerQAHistograms.cxx
  MpdTpcDigitizerTask.cxx
  MpdTpcHit.cxx
  MpdTpcHitProducer.cxx
  MpdTPCpid.cxx
  MpdTpcSector.cxx
  TaskHelpers.cxx
)

Set(HEADERS
  MpdParticleIdentification.h
  MpdTpc2dCluster.h
  MpdTpcClusterFinderMlem.h
  MpdTpcDigit.h
  MpdTpcDigitizerAZlt.h
  MpdTpcDigitizerQAHistograms.h
  MpdTpcDigitizerTask.h
  MpdTpcHit.h
  MpdTpcHitProducer.h
  MpdTPCpid.h
  MpdTpcSector.h
  TaskHelpers.h
)



Set(LINKDEF tpcLinkDef.h)
Set(LIBRARY_NAME tpc)
Set(DEPENDENCIES Base MpdField MpdMCStack Tof Kalman MpdGen tpcClusterHitFinder tpcDigitizer tpcGeometry tpcFairTpc)

GENERATE_LIBRARY()
