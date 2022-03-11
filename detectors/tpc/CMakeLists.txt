# Create a library called "libtpc" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_SOURCE_DIR}/core/mpdField
  ${CMAKE_SOURCE_DIR}/detectors/tof
  ${CMAKE_SOURCE_DIR}/generators
  ${CMAKE_SOURCE_DIR}/reconstruction/tracking/kalman
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
  MpdTpcClusterFinderAZ.cxx
  MpdTpcClusterFinderMlem.cxx
  MpdTpcClusterFinderQAHistograms.cxx
  MpdTpcClusterFinderTask.cxx
  MpdTpcDigit.cxx
  MpdTpcDigitizerAZ.cxx
  MpdTpcDigitizerAZlt.cxx
  MpdTpcDigitizerQAHistograms.cxx
  MpdTpcDigitizerTask.cxx
  MpdTpcEDepParams.cxx
  MpdTpcFoundHit.cxx
  MpdTpcHit.cxx
  MpdTpcHitProducer.cxx
  MpdTpcPeak.cxx
  MpdTPCpid.cxx
  MpdTpcSector.cxx
  MpdTpcSectorGeo.cxx
  TaskHelpers.cxx
  TpcContFact.cxx
  TpcDetector.cxx
  TpcGas.cxx
  TpcGeo.cxx
  TpcGeoPar.cxx
  TpcPoint.cxx
  TpcTimeBin.cxx
)

Set(HEADERS)
Set(LINKDEF tpcLinkDef.h)
Set(LIBRARY_NAME tpc)
Set(DEPENDENCIES Base Gdml MpdField MpdMCStack Tof Kalman MpdGen)

GENERATE_LIBRARY()
