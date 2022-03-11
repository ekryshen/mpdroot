# Create a library called "libTof" which includes the source files given in
# the array.
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_SOURCE_DIR}/kalman
  ${CMAKE_SOURCE_DIR}/simulation/mcStack
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
  MpdTof.cxx
  MpdTofBayesPid.cxx
  MpdTofContFact.cxx
  MpdTofGeo.cxx
  MpdTofGeoPar.cxx
  MpdTofGeoUtils.cxx
  MpdTofHit.cxx
  MpdTofHitProducer.cxx
  MpdTofHitProducerIdeal.cxx
  MpdTofHitProducerQA.cxx
  MpdTofMatching.cxx
  MpdTofMatchingData.cxx
  MpdTofMatchingQA.cxx
  MpdTOFpid.cxx
  MpdTofPoint.cxx
  MpdTofT0.cxx
)

set(HEADERS
  IntervalTree.h
  MpdTof.h
  MpdTofBayesPid.h
  MpdTofContFact.h
  MpdTofGeo.h
  MpdTofGeoPar.h
  MpdTofGeoUtils.h
  MpdTofHit.h
  MpdTofHitProducer.h
  MpdTofHitProducerIdeal.h
  MpdTofHitProducerQA.h
  MpdTofMatching.h
  MpdTofMatchingData.h
  MpdTofMatchingQA.h
  MpdTOFpid.h
  MpdTofPoint.h
  MpdTofT0.h
  MpdTofUtils.h
)

set(LINKDEF TofLinkDef.h)
set(LIBRARY_NAME Tof)
set(DEPENDENCIES MpdMCStack Kalman)

GENERATE_LIBRARY()
