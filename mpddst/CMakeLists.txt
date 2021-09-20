# Create a library called "libMpdDst" which includes the source files given in the array.
# The extension is already found. Any number of sources could be listed here.

add_subdirectory(mcDst)
add_subdirectory(MpdMiniEvent)

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_SOURCE_DIR}/mcstack
 ${CMAKE_SOURCE_DIR}/mpdbase
 ${CMAKE_SOURCE_DIR}/mpdfield
 ${CMAKE_SOURCE_DIR}/tpc
 ${CMAKE_SOURCE_DIR}/tof
 ${CMAKE_SOURCE_DIR}/emc
 ${CMAKE_SOURCE_DIR}/emc/emcKI
 ${CMAKE_SOURCE_DIR}/etof
 ${CMAKE_SOURCE_DIR}/kalman
 ${CMAKE_SOURCE_DIR}/lhetrack
 ${CMAKE_SOURCE_DIR}/mpdpid
 ${CMAKE_SOURCE_DIR}/zdc
 ${CMAKE_SOURCE_DIR}/mpddst
 ${CMAKE_SOURCE_DIR}/mpddst/MpdMiniEvent
)

Set(SYSTEM_INCLUDE_DIRECTORIES
 ${ROOT_INCLUDE_DIR}
)

include_directories(${INCLUDE_DIRECTORIES})
Include_Directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
 ${ROOT_LIBRARY_DIR}
 ${FAIRROOT_LIBRARY_DIR}
)

link_directories(${LINK_DIRECTORIES})

# List of source files
set(SRCS
 MpdFillDstTask.cxx
 MpdMiniDstFillTask.cxx
 MpdDstCompressTask.cxx
 MpdMiniEvent/MpdMiniBECalCluster.cxx
 MpdMiniEvent/MpdMiniBTofHit.cxx
 MpdMiniEvent/MpdMiniBTofPidTraits.cxx
 MpdMiniEvent/MpdMiniDstReader.cxx
 MpdMiniEvent/MpdMiniEvent.cxx
 MpdMiniEvent/MpdMiniFHCalHit.cxx
 MpdMiniEvent/MpdMiniHelix.cxx
 MpdMiniEvent/MpdMiniMcEvent.cxx
 MpdMiniEvent/MpdMiniMcTrack.cxx
 MpdMiniEvent/MpdMiniPhysicalHelix.cxx
 MpdMiniEvent/MpdMiniTrack.cxx
 MpdMiniEvent/MpdMiniTrackCovMatrix.cxx
)

Set(HEADERS)
Set(LINKDEF MpdDstLinkDef.h)
Set(LIBRARY_NAME MpdDst)
Set(DEPENDENCIES Core Base MpdBase LHETrack MpdPid MpdMiniEvent)

GENERATE_LIBRARY()
