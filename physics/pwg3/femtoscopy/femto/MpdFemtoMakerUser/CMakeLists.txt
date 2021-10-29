# Create a library called "MpdFemtoMakerUser"

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
 ${CMAKE_SOURCE_DIR}/physics/femto/MpdFemtoMakerUser
 ${CMAKE_SOURCE_DIR}/mpddst/mcDst
 ${CMAKE_SOURCE_DIR}/mpddst/MpdMiniEvent
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
MpdFemtoDummyTrackCut.cxx
MpdFemtoCorrFctn3DLCMSSym.cxx
MpdFemtoSmearPair.cxx
MpdFemtoModelQinvCorrFctn.cxx
MpdFemtoBPLCMS3DCorrFctnKt.cxx
MpdFemtoBasicEventCut.cxx
MpdFemtoBasicPairCut.cxx
MpdFemtoBasicTrackCut.cxx
MpdFemtoModelBPLCMS3DCorrFctnKt.cxx
MpdFemtoModelQinvCorrFctnKt.cxx
MpdFemtoQinvCorrFctnKt.cxx
MpdFemtoDeltaEtaDeltaPhiStarMinKt.cxx
MpdFemtoMcDstReader.cxx
MpdFemtoReactionPlaneAnalysis.cxx
MpdFemtoVertexAnalysis.cxx
MpdFemtoVertexMultAnalysis.cxx
MpdFemtoFxtEventCutMonitor.cxx
MpdFemtoFxtPairCutMonitor.cxx
MpdFemtoFxtTrackCutMonitor.cxx
MpdFemtoMiniDstReader.cxx
MpdFemtoAverageSeparation.cxx
        
MpdFemtoModelCorrFctnMomResolution.cxx        
)

Set(HEADERS)
Set(LINKDEF UserLinkDef.h)
Set(LIBRARY_NAME MpdFemtoMakerUser)
Set(DEPENDENCIES Core Base MpdBase MpdFemtoMaker MpdMcDst)

GENERATE_LIBRARY()
