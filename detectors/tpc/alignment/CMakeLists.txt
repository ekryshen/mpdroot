# Create a library called "libtpcAlignment"

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_CURRENT_LIST_DIR}
 ${CMAKE_CURRENT_LIST_DIR}/muons
 ${CMAKE_SOURCE_DIR}/detectors/tpc/geometry
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
  muons/TpcMissAlignment.cxx
)

Set(HEADERS
 muons/TpcMissAlignment.h
)

Set(LINKDEF tpcAlignmentLinkDef.h)
Set(LIBRARY_NAME tpcAlignment)
Set(DEPENDENCIES tpc MpdField MpdBase Minuit2)

GENERATE_LIBRARY()
