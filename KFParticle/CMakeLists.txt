cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)

#project(KFParticleLibrary)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/KFParticle/Vc/cmake" ${CMAKE_MODULE_PATH})



set(FIXTARGET FALSE CACHE BOOL "Compile for fix target geometry.")

#find_package(ROOT REQUIRED COMPONENTS Core Hist Matrix Physics EG Gpad Graf Graf3d RIO MathCore)
#If(DEFINED Vc_INCLUDE_DIR AND Vc_LIBRARIES)
#  Message("Vc found")
#Else()
#  find_package(Vc REQUIRED)
#EndIf()

Set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_SOURCE_DIR}/KFParticle
  ${CMAKE_SOURCE_DIR}/KFParticle/Vc
  ${CMAKE_SOURCE_DIR}/KFParticle/KFParticle
  ${CMAKE_SOURCE_DIR}/KFParticle/KFParticlePerformance

  ${CMAKE_SOURCE_DIR}/kalman
  ${CMAKE_SOURCE_DIR}/mcstack
  ${CMAKE_SOURCE_DIR}/mpdbase
  ${CMAKE_SOURCE_DIR}/mpddst
  ${CMAKE_SOURCE_DIR}/mpddst/MpdMiniEvent
)

Set(SYSTEM_INCLUDE_DIRECTORIES
 ${ROOT_INCLUDE_DIR} 
 ${CLHEP_INCLUDE_DIR}
 ${Boost_INCLUDE_DIRS}
)

set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)
Include_Directories( ${CMAKE_SOURCE_DIR} ${Vc_INCLUDE_DIR} ${INCLUDE_DIRECTORIES} ${SYSTEM_INCLUDE_DIRECTORIES})
set(CMAKE_INCLUDE_DIRECTORIES_BEFORE OFF)

set(LINK_DIRECTORIES
  ${ROOT_LIBRARY_DIR}
  ${FAIRROOT_LIBRARY_DIR}
)
 
link_directories( ${LINK_DIRECTORIES})

set (SRCS
  KFParticle/KFParticleTopoReconstructor.cxx
  KFParticle/KFVertex.cxx
  KFParticle/KFPTrack.cxx
  KFParticle/KFPTrackVector.cxx
  KFParticle/KFPVertex.cxx
  KFParticle/KFParticlePVReconstructor.cxx
  KFParticle/KFParticleDatabase.cxx
  KFParticle/KFParticleBase.cxx
  KFParticle/KFParticleBaseSIMD.cxx
  KFParticle/KFParticle.cxx
  KFParticle/KFParticleSIMD.cxx
  KFParticle/KFParticleFinder.cxx
  KFParticle/KFPEmcCluster.cxx

  KFParticlePerformance/KFMCVertex.cxx
  KFParticlePerformance/KFParticlePerformanceBase.cxx
  KFParticlePerformance/KFTopoPerformance.cxx
  KFParticlePerformance/KFMCParticle.cxx

  MpdKfParticleFinder.cxx
)

set (HEADERS
  KFParticle/KFParticleBase.h
  KFParticle/KFParticle.h
  KFParticle/KFVertex.h
  KFParticlePerformance/KFPartEfficiencies.h
  KFParticlePerformance/KFMCParticle.h

  MpdKfParticleFinder.h
)

Set(NODICT_HEADERS
  KFParticle/KFParticleTopoReconstructor.h
  KFParticle/KFParticlePVReconstructor.h
  KFParticle/KFPVertex.h
  KFParticle/KFPTrack.h
  KFParticle/KFParticleFinder.h
  KFParticle/KFParticleSIMD.h
  KFParticle/KFParticleBaseSIMD.h
  KFParticle/KFParticleDef.h
  KFParticle/KFPSimdAllocator.h
  KFParticle/KFPTrackVector.h
  KFParticle/KFPEmcCluster.h
  KFParticle/KFParticleField.h
  KFParticle/KFPTrackVector.h
  KFParticle/KFParticleDatabase.h
  KFParticlePerformance/KFTopoPerformance.h
  KFParticlePerformance/KFParticlePerformanceBase.h
  KFParticlePerformance/KFMCCounter.h
  KFParticlePerformance/KFPVEfficiencies.h
  KFParticlePerformance/KFMCVertex.h
  KFParticlePerformance/KFMCTrack.h
  KFParticlePerformance/KFPartMatch.h
)

Set(HEADERS)
Set(NODICT_HEADERS)
Set(LINKDEF KFLinkDef.h)
Set(LIBRARY_NAME MpdKfParticle)
Set(DEPENDENCIES Core Hist Matrix Physics EG Gpad Graf Graf3d RIO MathCore Vc Kalman MpdDst)

GENERATE_LIBRARY()
