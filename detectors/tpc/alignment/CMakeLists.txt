# Create a library called "libtpcAlignment"

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_CURRENT_LIST_DIR}
 ${CMAKE_CURRENT_LIST_DIR}/muons
 ${CMAKE_SOURCE_DIR}/detectors/tpc/
 ${CMAKE_SOURCE_DIR}/detectors/tpc/geometry
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
  #muons/TpcMissAlignment.cxx
  muons/TpcAlignment.cxx
  muons/Mille.cxx
  muons/DerLMpdTpc.cxx
  muons/DerHMpdTpc.cxx
  muons/RecoMiniMC.cxx
  #muons/TpcMiniMC.cxx
)

Set(HEADERS
  #muons/TpcMissAlignment.h
  muons/TpcAlignment.h
  muons/Mille.h
  muons/DerLMpdTpc.h
  muons/DerHMpdTpc.h
  muons/RecoMiniMC.h
  #muons/TpcMiniMC.h
)

# see alignment/muons/macro/ReadMe.md
if(ALIGNMENT_XI2)
    set (SRCS ${SRCS} muons/TpcMiniMC.cxx)
    set (HEADERS ${HEADERS} muons/TpcMiniMC.h)
else()
    set (SRCS ${SRCS} muons/TpcMissAlignment.cxx)
    set (HEADERS ${HEADERS} muons/TpcMissAlignment.h)
endif(ALIGNMENT_XI2)

Set(LINKDEF tpcAlignmentLinkDef.h)
Set(LIBRARY_NAME tpcAlignment)
Set(DEPENDENCIES tpc MpdField MpdBase Minuit2 tpcGeometry)

GENERATE_LIBRARY()
