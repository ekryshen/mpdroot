# Create a library called "libMpdDst" which includes the source files given in the array.
# The extension is already found. Any number of sources could be listed here.

add_subdirectory(mcDst)
add_subdirectory(MpdMiniEvent)

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_SOURCE_DIR}/mcstack
  ${CMAKE_SOURCE_DIR}/mpdbase
  ${CMAKE_SOURCE_DIR}/mpdfield
  ${CMAKE_SOURCE_DIR}/detectors/tpc
  ${CMAKE_SOURCE_DIR}/detectors/tof
  ${CMAKE_SOURCE_DIR}/detectors/emc/emcKI
  ${CMAKE_SOURCE_DIR}/kalman
  ${CMAKE_SOURCE_DIR}/lhetrack
  ${CMAKE_SOURCE_DIR}/mpdpid
  ${CMAKE_SOURCE_DIR}/detectors/zdc
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
)

Set(HEADERS)
Set(LINKDEF MpdDstLinkDef.h)
Set(LIBRARY_NAME MpdDst)
Set(DEPENDENCIES MpdMiniEvent)

GENERATE_LIBRARY()
