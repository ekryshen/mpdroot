# Create a library called "MpdFemtoMaker"

set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_SOURCE_DIR}/physics/femto/MpdFemtoMaker
 )

set(SYSTEM_INCLUDE_DIRECTORIES
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

Set(HEADERS
  math_constants.h
  MpdFemtoAnalysisCollection.h
  MpdFemtoBaseAnalysis.h
  MpdFemtoBaseEventCut.h
  MpdFemtoBaseEventWriter.h
  MpdFemtoBaseKinkCut.h
  MpdFemtoBaseLikeSignCorrFctn.h
  MpdFemtoBasePairCut.h
  MpdFemtoBaseParticleCut.h
  MpdFemtoBaseTrackCut.h
  MpdFemtoBaseV0Cut.h
  MpdFemtoBaseXiCut.h
  MpdFemtoBaseCorrFctn.h
  MpdFemtoCorrFctnCollection.h
  MpdFemtoCutMonitorCollection.h
  MpdFemtoEnumeration.h
  MpdFemtoKink.h
  MpdFemtoLikeSignAnalysis.h
  MpdFemtoPair.h
  MpdFemtoEvent.h
  MpdFemtoEventCutCollection.h
  MpdFemtoBaseCutMonitor.h
  MpdFemtoModelHiddenInfo.h
  MpdFemtoEventWriterCollection.h
  MpdFemtoParticle.h
  MpdFemtoPhysicalHelix.h
  MpdFemtoHiddenInfo.h
  MpdFemtoBaseEventReader.h
  MpdFemtoKinkCollection.h
  MpdFemtoManager.h
  MpdFemtoModelManager.h
  MpdFemtoAnalysis.h
  MpdFemtoCutMonitorHandler.h
  MpdFemtoHelix.h
  MpdFemtoBaseModelFreezeOutGenerator.h
  MpdFemtoBaseModelWeightGenerator.h
  MpdFemtoModelWeightGeneratorLednicky.h
  MpdFemtoModelGausLCMSFreezeOutGenerator.h
  MpdFemtoPicoEventCollectionVectorHideAway.h
  MpdFemtoMultiTrackCut.h
  MpdFemtoMaker.h
  MpdFemtoPairCutCollection.h
  MpdFemtoParticleCollection.h
  MpdFemtoPicoEvent.h
  MpdFemtoTriplet.h
  MpdFemtoPicoEventCollection.h
  MpdFemtoPicoEventCollectionVector.h
  MpdFemtoString.h
  MpdFemtoTrack.h
  MpdFemtoTrackCollection.h
  MpdFemtoTrackCutCollection.h
  MpdFemtoTypes.h
  MpdFemtoV0.h
  MpdFemtoV0Collection.h
  MpdFemtoXiCollection.h
  MpdFemtoXi.h
  TpcLocalTransform.h
  phys_constants.h
  PhysicalConstants.h
  SystemOfUnits.h
)

set(LINKDEF BaseLinkDef.h)
set(LIBRARY_NAME MpdFemtoMaker)
set(DEPENDENCIES Base)

GENERATE_LIBRARY()
