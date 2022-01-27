
If(DEFINED Vc_INCLUDE_DIR AND Vc_LIBRARIES)
  Message("Vc found")
Else()
  find_package(Vc REQUIRED)
EndIf()


message(WARNING "===== IN KALMAN")
Set(INCLUDE_DIRECTORIES
  ${CMAKE_SOURCE_DIR}/external/kfparticle
  ${CMAKE_SOURCE_DIR}/external/kfparticle/KFParticle
  ${CMAKE_SOURCE_DIR}/external/kfparticle/KFParticlePerformance
)

Include_Directories( ${CMAKE_SOURCE_DIR} ${Vc_INCLUDE_DIR} ${INCLUDE_DIRECTORIES} )

set(LINK_DIRECTORIES
  ${ROOT_LIBRARY_DIR}
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
  KFParticleTest/KFParticleTest.cxx
)
#  ROOT_GENERATE_DICTIONARY(G__KFParticle ${HEADERS} LINKDEF KFLinkDef.h OPTIONS "-DDO_TPCCATRACKER_EFF_PERFORMANCE" "-DHomogeneousField" "-DUSE_TIMERS")
#  add_library(KFParticle SHARED ${SOURCES} G__KFParticle.cxx)
#  target_link_libraries(KFParticle ${ROOT_LIBRARIES} ${Vc_LIBRARIES} )
# add_target_property(KFParticle COMPILE_FLAGS "-DDO_TPCCATRACKER_EFF_PERFORMANCE -DHomogeneousField -DUSE_TIMERS")
  
  
#add_definitions(-DDO_TPCCATRACKER_EFF_PERFORMANCE -DHomogeneousField -DUSE_TIMERS)
#add_compile_options("-DDO_TPCCATRACKER_EFF_PERFORMANCE -DHomogeneousField -DUSE_TIMERS")

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

set (HEADERS
  KFParticle/KFParticleBase.h
  KFParticle/KFParticle.h
  KFParticle/KFVertex.h
  KFParticlePerformance/KFPartEfficiencies.h
  KFParticlePerformance/KFMCParticle.h
  KFParticleTest/KFParticleTest.h
)

Set(LINKDEF KFLinkDef.h)
Set(LIBRARY_NAME KFParticle)
Set(DEPENDENCIES Core Hist Matrix Physics EG Gpad Graf Graf3d RIO MathCore  Vc)
set(DEFINITIONS -DDO_TPCCATRACKER_EFF_PERFORMANCE -DHomogeneousField -DUSE_TIMERS)
#ROOT_GENERATE_DICTIONARY(G__KFParticle ${HEADERS} LINKDEF KFLinkDef.h OPTIONS "-DDO_TPCCATRACKER_EFF_PERFORMANCE" "-DHomogeneousField" "-DUSE_TIMERS")
#  add_library(KFParticle SHARED ${SOURCES} G__KFParticle.cxx)
#  target_link_libraries(KFParticle ${ROOT_LIBRARIES} ${Vc_LIBRARIES} )
  
#  add_target_property(KFParticle COMPILE_FLAGS "-DDO_TPCCATRACKER_EFF_PERFORMANCE -DHomogeneousField -DUSE_TIMERS")
#install(FILES ${HEADERS} ${NODICT_HEADERS} DESTINATION include)
#install(FILES ${CMAKE_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}KFParticle_rdict.pcm DESTINATION lib OPTIONAL)
#install(FILES ${CMAKE_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}KFParticle.rootmap DESTINATION lib)
#install(TARGETS KFParticle DESTINATION lib)
GENERATE_LIBRARY()
