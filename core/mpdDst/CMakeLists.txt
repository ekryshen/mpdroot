# Create a library called "libMpdDst" which includes the source files given in the array.
# The extension is already found. Any number of sources could be listed here.
add_subdirectory(MpdMiniEvent)

set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_CURRENT_LIST_DIR}/MpdMiniEvent
  ${CMAKE_SOURCE_DIR}/core/mpdBase
  ${CMAKE_SOURCE_DIR}/core/mpdField
  ${CMAKE_SOURCE_DIR}/core/mpdPid
  ${CMAKE_SOURCE_DIR}/detectors/emc/emcKI
  ${CMAKE_SOURCE_DIR}/detectors/multi
  ${CMAKE_SOURCE_DIR}/detectors/tpc
  ${CMAKE_SOURCE_DIR}/detectors/tpc/digitizer/AZlt
  ${CMAKE_SOURCE_DIR}/detectors/tpc/clusterHitFinder/mlem
  ${CMAKE_SOURCE_DIR}/detectors/tpc/clusterHitFinder/fast
  ${CMAKE_SOURCE_DIR}/detectors/tpc/pid
  ${CMAKE_SOURCE_DIR}/detectors/tof
  ${CMAKE_SOURCE_DIR}/detectors/zdc
  ${CMAKE_SOURCE_DIR}/reconstruction/tracking/kalman
  ${CMAKE_SOURCE_DIR}/reconstruction/tracking/lheTrack
  ${CMAKE_SOURCE_DIR}/simulation/mcStack
  ${CMAKE_SOURCE_DIR}/tools/tdd/QA
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
Set(DEPENDENCIES MpdMiniEvent multi)

GENERATE_LIBRARY()
