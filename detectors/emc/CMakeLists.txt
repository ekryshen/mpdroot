# Create a library called "libEmc" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_CURRENT_LIST_DIR}
 ${CMAKE_CURRENT_LIST_DIR}/hitmaker
 ${CMAKE_CURRENT_LIST_DIR}/emcKI
 ${CMAKE_SOURCE_DIR}/clustering
 ${CMAKE_SOURCE_DIR}/mcstack
 ${CMAKE_SOURCE_DIR}/detectors/tpc
 ${CMAKE_SOURCE_DIR}/kalman
 ${CMAKE_SOURCE_DIR}/mpdbase
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
MpdEmcContFact.cxx
MpdEmc.cxx
MpdEmcGeo.cxx
MpdEmcGeoPar.cxx
MpdEmcGeoParams.cxx
MpdEmcPoint.cxx
MpdEmcHit.cxx
MpdEmcCluster.cxx
MpdEmcHitProducer.cxx
MpdEmcHitCreation.cxx
MpdEmcClusterCreation.cxx

MpdEmcDigit.cxx
MpdEmcDigitizer.cxx
MpdEmcClusterFinderAZ.cxx
MpdEmcMatching.cxx
MpdEmcMatch.cxx
MpdEmcCalibParams.cxx

hitmaker/MpdEmcHitA.cxx
hitmaker/MpdEmcHitProducerA.cxx
hitmaker/MpdEmcHitCont.cxx
hitmaker/MpdEmcGeoParWrapper.cxx

clustering/EmcClManager.cxx
clustering/EmcClusterElement.cxx
clustering/EmcCluster.cxx
clustering/EmcClusterCollection.cxx
clustering/EmcClusterInfo.cxx
clustering/EmcCLinfo.cxx

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
)

Set(HEADERS)
Set(LINKDEF MpdEmcLinkDef.h)
Set(LIBRARY_NAME Emc)
Set(DEPENDENCIES)

GENERATE_LIBRARY()
