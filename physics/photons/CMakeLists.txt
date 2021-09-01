set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_SOURCE_DIR}/mpdbase
 ${CMAKE_SOURCE_DIR}/mpdpid
 ${CMAKE_SOURCE_DIR}/kalman
 ${CMAKE_SOURCE_DIR}/lhetrack
 ${CMAKE_SOURCE_DIR}/emc/emcKI
 ${CMAKE_SOURCE_DIR}/mcstack
 ${CMAKE_SOURCE_DIR}/physics
 ${CMAKE_SOURCE_DIR}/physics/photons
)

Set(SYSTEM_INCLUDE_DIRECTORIES
 ${ROOT_INCLUDE_DIR}
 ${Boost_INCLUDE_DIRS}
)

include_directories(${INCLUDE_DIRECTORIES})
Include_Directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
 ${ROOT_LIBRARY_DIR}
 ${FAIRROOT_LIBRARY_DIR}
) 

link_directories(${LINK_DIRECTORIES})

set(SRCS
  MpdPhoton.h
  MpdConvPi0.cxx  
  MpdPhotonAnalysisParams.cxx 
)

Set(HEADERS)


Set(LINKDEF PhotonsLinkDef.h)
Set(LIBRARY_NAME MpdPhotons)
Set(DEPENDENCIES MpdBase MpdMCStack Kalman MpdPid LHETrack Physics)

GENERATE_LIBRARY()
