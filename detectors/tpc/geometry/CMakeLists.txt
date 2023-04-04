# Create a library called "libtpcAlignment"

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_CURRENT_LIST_DIR}
 ${CMAKE_SOURCE_DIR}/detectors/tpc
 ${CMAKE_SOURCE_DIR}/detectors/tpc/fairTpc
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
  TpcSectorGeoAlignmentVAK.cxx
  TpcSectorGeoAZ.cxx
)

Set(HEADERS
TpcSectorGeoAlignmentVAK.h
TpcSectorGeoAZ.h
)

Set(LINKDEF tpcGeometryLinkDef.h)
Set(LIBRARY_NAME tpcGeometry)
Set(DEPENDENCIES tpc MpdField MpdBase tpcFairTpc)

GENERATE_LIBRARY()
