# Create a library called "libtpcAlignment"

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_CURRENT_LIST_DIR}
 ${CMAKE_SOURCE_DIR}/detectors/tpc
 ${CMAKE_SOURCE_DIR}/detectors/tpc/fairTpc
 ${CMAKE_SOURCE_DIR}/detectors/tpc/geometry
 ${CMAKE_SOURCE_DIR}/simulation/mcStack
 ${CMAKE_SOURCE_DIR}/core/mpdField
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
 AZlt/MpdTpcDigit.cxx
 AZlt/MpdTpcDigitizerAZlt.cxx
 AZlt/MpdTpcDigitizerQAHistograms.cxx
 AZlt/MpdTpcDigitizerTask.cxx
 AZlt/MpdTpcSector.cxx
 AZlt/TaskHelpers.cxx
)

Set(HEADERS
 AbstractTpcDigit.h
 AZlt/MpdTpcDigit.h
 AZlt/MpdTpcDigitizerAZlt.h
 AZlt/MpdTpcDigitizerQAHistograms.h
 AZlt/MpdTpcDigitizerTask.h
 AZlt/MpdTpcSector.h
 AZlt/TaskHelpers.h
)

Set(LINKDEF tpcDigitizerLinkDef.h)
Set(LIBRARY_NAME tpcDigitizer)
Set(DEPENDENCIES MpdMCStack MpdField)

GENERATE_LIBRARY()
