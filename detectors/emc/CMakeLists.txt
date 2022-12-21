# Create a library called "libEmc" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_CURRENT_LIST_DIR}/emcClustering
  ${CMAKE_CURRENT_LIST_DIR}/emcClustering/clustering
  ${CMAKE_CURRENT_LIST_DIR}/emcKI
  ${CMAKE_CURRENT_LIST_DIR}/hitMaker
  ${CMAKE_SOURCE_DIR}/core/mpdBase
  ${CMAKE_SOURCE_DIR}/detectors/tpc
  ${CMAKE_SOURCE_DIR}/detectors/tpc/clusterHitFinder
  ${CMAKE_SOURCE_DIR}/detectors/tpc/digitizer
  ${CMAKE_SOURCE_DIR}/reconstruction/tracking/kalman
  ${CMAKE_SOURCE_DIR}/simulation/mcStack
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
  MpdEmc.cxx
  MpdEmcCalibParams.cxx
  MpdEmcCluster.cxx
  MpdEmcClusterCreation.cxx
  MpdEmcClusterFinderAZ.cxx
  MpdEmcContFact.cxx
  MpdEmcDigit.cxx
  MpdEmcDigitizer.cxx
  MpdEmcGeo.cxx
  MpdEmcGeoPar.cxx
  MpdEmcGeoParams.cxx
  MpdEmcHit.cxx
  MpdEmcHitCreation.cxx
  MpdEmcHitProducer.cxx
  MpdEmcMatch.cxx
  MpdEmcMatching.cxx
  MpdEmcPoint.cxx

  emcClustering/clustering/ClusterElement.cxx
  emcClustering/clustering/Cluster.cxx
  emcClustering/clustering/ClusterCollection.cxx
  emcClustering/clustering/ClManager.cxx

  emcClustering/EmcCLinfo.cxx
  emcClustering/EmcClManager.cxx
  emcClustering/EmcCluster.cxx
  emcClustering/EmcClusterCollection.cxx
  emcClustering/EmcClusterElement.cxx
  emcClustering/EmcClusterInfo.cxx

  emcKI/MpdEmcClusterizerKI.cxx
  emcKI/MpdEmcClusterKI.cxx
  emcKI/MpdEmcDigitizerKI.cxx
  emcKI/MpdEmcDigitKI.cxx
  emcKI/MpdEmcGeoUtils.cxx
  emcKI/MpdEmcKI.cxx
  emcKI/MpdEmcMatchingKI.cxx
  emcKI/MpdEmcPointKI.cxx
  emcKI/MpdEmcSimParams.cxx
  emcKI/MpdEmcTrackExtrap.cxx

  hitMaker/MpdEmcGeoParWrapper.cxx
  hitMaker/MpdEmcHitA.cxx
  hitMaker/MpdEmcHitCont.cxx
  hitMaker/MpdEmcHitProducerA.cxx
)

set(HEADERS
#  emc.h
  MpdEmc.h
  MpdEmcCalibParams.h
  MpdEmcCluster.h
  MpdEmcClusterCreation.h
  MpdEmcClusterFinderAZ.h
  MpdEmcContFact.h
  MpdEmcDigit.h
  MpdEmcDigitizer.h
  MpdEmcGeo.h
  MpdEmcGeoPar.h
  MpdEmcGeoParams.h
  MpdEmcHit.h
  MpdEmcHitCreation.h
  MpdEmcHitProducer.h
  MpdEmcMatch.h
  MpdEmcMatching.h
  MpdEmcPoint.h

  emcClustering/clustering/ClusterElement.h
  emcClustering/clustering/Cluster.h
  emcClustering/clustering/ClusterCollection.h
  emcClustering/clustering/ClManager.h

  emcClustering/EmcCLinfo.h
  emcClustering/EmcClManager.h
  emcClustering/EmcCluster.h
  emcClustering/EmcClusterCollection.h
  emcClustering/EmcClusterElement.h
  emcClustering/EmcClusterInfo.h

  emcKI/MpdEmcClusterizerKI.h
  emcKI/MpdEmcClusterKI.h
  emcKI/MpdEmcDigitizerKI.h
  emcKI/MpdEmcDigitKI.h
  emcKI/MpdEmcGeoUtils.h
  emcKI/MpdEmcKI.h
  emcKI/MpdEmcMatchingKI.h
  emcKI/MpdEmcPointKI.h
  emcKI/MpdEmcSimParams.h
  emcKI/MpdEmcTrackExtrap.h

  hitMaker/MpdEmcGeoParWrapper.h
  hitMaker/MpdEmcHitA.h
  hitMaker/MpdEmcHitCont.h
  hitMaker/MpdEmcHitProducerA.h
)

Set(LINKDEF MpdEmcLinkDef.h)
Set(LIBRARY_NAME Emc)
Set(DEPENDENCIES MpdMCStack tpc Kalman MpdBase)

GENERATE_LIBRARY()
