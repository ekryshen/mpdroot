# Create a library called "libtpc" which is an API for MPD TPC detector

set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_SOURCE_DIR}/tools/tdd/QA
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
 AbstractTpcClusterHitFinder.cxx
 BaseTpcSectorGeo.cxx
)

Set(HEADERS
 AbstractTpc2dCluster.h
 AbstractTpcClusterHitFinder.h
 AbstractTpcHit.h
 AbstractTpcDigit.h
 BaseTpcSectorGeo.h
)


Set(LINKDEF tpcLinkDef.h)
Set(LIBRARY_NAME tpc)
Set(DEPENDENCIES QA)

GENERATE_LIBRARY()
