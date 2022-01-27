MACRO(InstallKFPackage)
include(ExternalProject)
set(KFPARTICLE_LIBNAME "${CMAKE_SHARED_LIBRARY_PREFIX}KFParticle${CMAKE_SHARED_LIBRARY_SUFFIX}")

set(KFPARTICLE_SRC_URL "https://github.com/Ivan-Kisel/KFParticle.git")
set(KFPARTICLE_SRC "${CMAKE_SOURCE_DIR}/external/kfparticle")
set(KFPARTICLE_DESTDIR "${CMAKE_BINARY_DIR}/external/kfparticle")
set(KFPARTICLE_VER "v1.1")
# GIT_TAG is a hash for KFParticle tag cbm/v1.1-1 


#

if(NOT EXISTS "${KFPARTICLE_SRC}/CMakeLists.txt")
    execute_process(COMMAND git clone -b ${KFPARTICLE_VER} ${KFPARTICLE_SRC_URL} ${KFPARTICLE_SRC})
    execute_process(COMMAND git reset --hard  ${KFPARTICLE_VER} WORKING_DIRECTORY  ${KFPARTICLE_SRC})
endif()
execute_process(COMMAND git fetch WORKING_DIRECTORY  ${KFPARTICLE_SRC})
 #   execute_process(COMMAND git reset --hard  ${KFPARTICLE_VER} WORKING_DIRECTORY  ${CMAKE_SOURCE_DIR}/external/kfparticle)
        #add_subdirectory (${CMAKE_SOURCE_DIR}/external/kfpackage)


If(ProjectUpdated)
  File(REMOVE_RECURSE ${KFPARTICLE_SRC})
  Message("KFParticle source directory was changed so build directory was deleted")  
EndIf()
#add_subdirectory(${KFPARTICLE_SRC})
add_definitions(-DDO_TPCCATRACKER_EFF_PERFORMANCE -DHomogeneousField -DUSE_TIMERS)

file(COPY ${CMAKE_SOURCE_DIR}/cmake/external_patches/kfpackage/CMakeLists.txt DESTINATION ${KFPARTICLE_SRC})
add_subdirectory(external/kfparticle)


#set(KFParticle_LIB_DIR ${CMAKE_BINARY_DIR}/lib ${CMAKE_SOURCE_DIR}/build/lib)
#set(KFParticle_LIBRARIES KFParticle ${CMAKE_SOURCE_DIR}/build)
#set(KFParticle_INCLUDE_DIR "${CMAKE_BINARY_DIR}/include" ${CMAKE_SOURCE_DIR}/build/include)
set(KFParticle_FOUND TRUE ${CMAKE_SOURCE_DIR}/external/kfparticle)

Install(FILES ${CMAKE_BINARY_DIR}/lib/${KFPARTICLE_LIBNAME} 
              ${CMAKE_BINARY_DIR}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}KFParticle.rootmap
              ${CMAKE_BINARY_DIR}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}KFParticle_rdict.pcm
        DESTINATION lib
       )

Install(FILES ${CMAKE_BINARY_DIR}/include/KFParticleBase.h 
              ${CMAKE_BINARY_DIR}/include/KFParticle.h
              ${CMAKE_BINARY_DIR}/include/KFVertex.h
        DESTINATION include/KFParticle
       )
Install(FILES ${CMAKE_BINARY_DIR}/include/KFMCParticle.h
              ${CMAKE_BINARY_DIR}/include/KFPartEfficiencies.h
        DESTINATION include/KFParticlePerformance
       )
Install(FILES ${CMAKE_BINARY_DIR}/include/KFParticleTest.h
        DESTINATION include/KFParticleTest
       )
Install(FILES ${CMAKE_BINARY_DIR}/include/KFPVertex.h
              ${CMAKE_BINARY_DIR}/include/KFMCCounter.h
        DESTINATION include
       )
ENDMACRO(InstallKFPackage)
