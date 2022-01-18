

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_SOURCE_DIR}/mpddst
  ${CMAKE_SOURCE_DIR}/mpdbase
  ${CMAKE_SOURCE_DIR}/mpddst/mcDst
  ${CMAKE_SOURCE_DIR}/mpddst/MpdMiniEvent
  ${CMAKE_SOURCE_DIR}/physics/common/v0
  ${CMAKE_SOURCE_DIR}/external/nicafemto/features
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
	MpdV0DaughterCut.cxx
	MpdV0DaughterCutBasic.cxx
	MpdV0FinderBasic.cxx
    MpdV0FinderHelix.cxx
    MpdV0CandidateCut.cxx
	MpdV0Particle.cxx
	MpdV0Track.cxx
	MpdV0Namespace.cxx
)

Set(LINKDEF MpdV0CommonLinkDef.h)
Set(LIBRARY_NAME MpdV0Common)
Set(DEPENDENCIES Core Base MpdDst MpdBase NicaFeatures)

GENERATE_LIBRARY()