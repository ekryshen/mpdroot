add_subdirectory(femto/fsiTools)



add_subdirectory(femto/MpdFemtoMaker)



add_subdirectory(femto/MpdFemtoMakerUser)



add_subdirectory(photons)



add_subdirectory(pairKK)



add_subdirectory(pairPK)



add_subdirectory(pairPiK)



add_subdirectory(pairPiPi)



add_subdirectory(pairPiLambda)



add_subdirectory(pairPiKs)



add_subdirectory(evCentrality)



add_subdirectory(evPlane)



add_subdirectory(evPID)



add_subdirectory(globalPolarization)

add_subdirectory(hyperons)





# Create a library called "MpdPhysics" which includes the source files given in the array.



# The extension is already found. Any number of sources could be listed here.



set(INCLUDE_DIRECTORIES



 ${BASE_INCLUDE_DIRECTORIES}



  ${CMAKE_SOURCE_DIR}/core/mpdBase



  ${CMAKE_SOURCE_DIR}/physics



  ${CMAKE_SOURCE_DIR}/reconstruction/tracking/kalman



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



 MpdAnalysisEvent.cxx



 MpdAnalysisEventPlane.cxx



 MpdAnalysisTask.cxx



 MpdAnalysisManager.cxx



)







Set(HEADERS)



Set(LINKDEF MpdPhysicsLinkDef.h)



Set(LIBRARY_NAME MpdPhysics)



Set(DEPENDENCIES Base MpdBase)







GENERATE_LIBRARY()



