# Create a library called "libLHETrack" 

set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_SOURCE_DIR}/core/mpdBase
  ${CMAKE_SOURCE_DIR}/core/mpdField 
  ${CMAKE_SOURCE_DIR}/detectors/sts
  ${CMAKE_SOURCE_DIR}/detectors/tpc
  ${CMAKE_SOURCE_DIR}/detectors/tpc/geometry
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
  MpdEctKalmanTrack.cxx
  MpdItsKalmanTrack.cxx
  MpdItsToTpcMatching.cxx
  MpdKfPrimaryVertexFinder.cxx
  MpdMotherFitterPart.cxx
  MpdParticle.cxx
  MpdTpcDedxTask.cxx
  MpdTpcKalmanFilter.cxx
  MpdTrackFinderIts5spd.cxx
  MpdVector.cxx
  MpdVectorFinder.cxx
  MpdVertexZfinder.cxx
)

# List of header files
set(HEADERS)

Set(LINKDEF LHETrackLinkDef.h)
Set(LIBRARY_NAME LHETrack)
Set(DEPENDENCIES Kalman MpdBase)

GENERATE_LIBRARY()
