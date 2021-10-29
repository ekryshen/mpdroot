# Create a library called "MpdFemtoMaker"

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_SOURCE_DIR}/mcstack
 ${CMAKE_SOURCE_DIR}/passive
 ${CMAKE_SOURCE_DIR}/mpdbase
 ${CMAKE_SOURCE_DIR}/mpddst
 ${CMAKE_SOURCE_DIR}/mpdfield
 ${CMAKE_SOURCE_DIR}/tpc
 ${CMAKE_SOURCE_DIR}/lhetrack
 ${CMAKE_SOURCE_DIR}/kalman
 ${CMAKE_SOURCE_DIR}/physics/femto/MpdFemtoMaker
 )

Set(SYSTEM_INCLUDE_DIRECTORIES
 ${ROOT_INCLUDE_DIR}
 ${GSL_INCLUDE_DIRS}
 ${CLHEP_INCLUDE_DIR}
 ${Boost_INCLUDE_DIRS}
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
MpdFemtoMultiTrackCut.cxx
TpcLocalTransform.cxx
MpdFemtoMaker.cxx
MpdFemtoModelWeightGeneratorLednicky.cxx
MpdFemtoManager.cxx
MpdFemtoModelGausLCMSFreezeOutGenerator.cxx
MpdFemtoXi.cxx
MpdFemtoModelManager.cxx
MpdFemtoHelix.cxx
MpdFemtoCutMonitorHandler.cxx
MpdFemtoBaseCutMonitor.cxx
MpdFemtoModelHiddenInfo.cxx
MpdFemtoBaseModelWeightGenerator.cxx
MpdFemtoPicoEvent.cxx
MpdFemtoBaseModelFreezeOutGenerator.cxx
MpdFemtoTriplet.cxx
MpdFemtoPicoEventCollectionVectorHideAway.cxx
MpdFemtoAnalysis.cxx
MpdFemtoBaseCorrFctn.cxx
MpdFemtoBaseEventReader.cxx
MpdFemtoEvent.cxx
MpdFemtoKink.cxx
MpdFemtoLikeSignAnalysis.cxx
MpdFemtoPair.cxx
MpdFemtoParticle.cxx
MpdFemtoPhysicalHelix.cxx
MpdFemtoTrack.cxx
MpdFemtoV0.cxx
)

Set(HEADERS)
Set(LINKDEF BaseLinkDef.h)
Set(LIBRARY_NAME MpdFemtoMaker)
Set(DEPENDENCIES Core Base MpdBase MpdFsiTools)

GENERATE_LIBRARY()
