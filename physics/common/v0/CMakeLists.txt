

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_SOURCE_DIR}/mpddst
  ${CMAKE_SOURCE_DIR}/mpdbase
  ${CMAKE_SOURCE_DIR}/kalman
  ${CMAKE_SOURCE_DIR}/mcstack/
  ${CMAKE_SOURCE_DIR}/mpddst/mcDst
  ${CMAKE_SOURCE_DIR}/mpddst/MpdMiniEvent
  ${CMAKE_SOURCE_DIR}/physics/common/v0
  ${CMAKE_SOURCE_DIR}/physics/common/v0/cuts
    ${CMAKE_SOURCE_DIR}/physics/common/v0/monitors
  ${CMAKE_SOURCE_DIR}/external/nicafemto/features
  ${CMAKE_SOURCE_DIR}/external/kfparticle/KFParticle
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
	cuts/MpdV0DaughterCut.cxx
	cuts/MpdV0DaughterCutBasic.cxx
	cuts/MpdV0CandidateCut.cxx
	cuts/MpdV0CandidateCutBasic.cxx
	cuts/MpdV0CandidateCutKalman.cxx
	cuts/MpdV0FinderCut.cxx
	
	monitors/MpdV0CandidateMonitor.cxx
	monitors/MpdV0DaughterMonitor.cxx
	monitors/MpdV0Monitor.cxx
	monitors/MpdV0StandardCandidateMonitor.cxx
	monitors/MpdV0StandardDaughterMonitor.cxx
	
	MpdV0Finder.cxx
    MpdV0FinderHelix.cxx
    MpdV0FinderKFPackage.cxx
    MpdV0Matcher.cxx
    
	MpdV0Particle.cxx
	MpdV0Track.cxx
	MpdV0Namespace.cxx
	MpdSimpleLinks.cxx
	MpdDCAOnTheFly.cxx
)

Set(LINKDEF MpdV0CommonLinkDef.h)
Set(LIBRARY_NAME MpdV0Common)
Set(DEPENDENCIES Core Base MpdDst MpdMiniEvent KFParticle Kalman MpdMCStack)

GENERATE_LIBRARY()