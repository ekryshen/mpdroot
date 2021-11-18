# Create a library called "libKalman" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_SOURCE_DIR}/mpdbase
 ${CMAKE_SOURCE_DIR}/mpdfield
 ${CMAKE_SOURCE_DIR}/kalman
)

Set(SYSTEM_INCLUDE_DIRECTORIES
 ${ROOT_INCLUDE_DIR}
)

include_directories(${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
 ${ROOT_LIBRARY_DIR}
 ${FAIRROOT_LIBRARY_DIR}
)
 
link_directories( ${LINK_DIRECTORIES})

# List of source files
set(SRCS
 MpdKalmanFilter.cxx
 MpdKalmanGeoScheme.cxx
 MpdKalmanHit.cxx
 MpdKalmanTrack.cxx
 MpdTpcKalmanTrack.cxx
)

Set(HEADERS)
Set(LINKDEF KalmanLinkDef.h)
Set(LIBRARY_NAME Kalman)
Set(DEPENDENCIES MpdField)

GENERATE_LIBRARY()
