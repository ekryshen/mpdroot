# Create a library called "libtpc" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_SOURCE_DIR}/mpdfield
 ${CMAKE_SOURCE_DIR}/mcstack
 ${CMAKE_SOURCE_DIR}/passive
 ${CMAKE_SOURCE_DIR}/tpc
 ${CMAKE_SOURCE_DIR}/tof
 ${CMAKE_SOURCE_DIR}/lhetrack
 ${CMAKE_SOURCE_DIR}/kalman
 ${CMAKE_SOURCE_DIR}/floor
)

Set(SYSTEM_INCLUDE_DIRECTORIES
 ${ROOT_INCLUDE_DIR} 
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
MpdFloor.cxx
MpdFloorGeo.cxx
MpdFloorGeoPar.cxx
MpdFloorPoint.cxx
MpdFloorContFact.cxx
)

Set(HEADERS)
Set(LINKDEF MpdFloorLinkDef.h)
Set(LIBRARY_NAME Floor)
Set(DEPENDENCIES Core Base)

GENERATE_LIBRARY()
