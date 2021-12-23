# Create a library called "libTof" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_SOURCE_DIR}/mcstack
  ${CMAKE_SOURCE_DIR}/tof
  ${CMAKE_SOURCE_DIR}/kalman
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
MpdTofGeo.cxx 
MpdTofGeoPar.cxx 
MpdTof.cxx 
MpdTofContFact.cxx 
MpdTofHit.cxx 
MpdTofHitProducer.cxx 
MpdTofHitProducerIdeal.cxx 
MpdTofPoint.cxx
MpdTofGeoUtils.cxx
MpdTofMatching.cxx
MpdTofMatchingData.cxx
MpdTofMatchingQA.cxx
MpdTofHitProducerQA.cxx
MpdTofT0.cxx
MpdTofBayesPid.cxx
MpdTOFpid.cxx
)

# List of header files
set(HEADERS 
MpdTofUtils.h
IntervalTree.h
MpdTofGeo.h
MpdTofGeoPar.h
MpdTof.h
MpdTofContFact.h
MpdTofHit.h
MpdTofHitProducer.h
MpdTofHitProducerIdeal.h
MpdTofPoint.h
MpdTofGeoUtils.h
MpdTofMatching.h
MpdTofMatchingData.h
MpdTofMatchingQA.h
MpdTofHitProducerQA.h
MpdTofT0.h
MpdTofBayesPid.h
MpdTOFpid.h
)

Set(LINKDEF TofLinkDef.h)

Set(LIBRARY_NAME Tof)
Set(DEPENDENCIES)

GENERATE_LIBRARY()
