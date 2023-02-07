# Create a library called "libtpcFairTpc"
# this library is old and will be replaced

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_CURRENT_LIST_DIR}
 ${CMAKE_SOURCE_DIR}/detectors/tpc
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
 TpcDetector.cxx
 TpcGas.cxx 
 TpcGeo.cxx 
 TpcGeoPar.cxx 
 TpcPoint.cxx 
)

Set(HEADERS
 TpcDetector.h
 TpcGas.h 
 TpcGeo.h 
 TpcGeoPar.h 
 TpcPoint.h 
)

Set(LINKDEF tpcFairTpcLinkDef.h)
Set(LIBRARY_NAME tpcFairTpc)
Set(DEPENDENCIES Gdml MpdMCStack)

GENERATE_LIBRARY()
