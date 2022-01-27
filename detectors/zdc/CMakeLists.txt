# Create a library called "libZdc" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_CURRENT_LIST_DIR}
 ${CMAKE_SOURCE_DIR}/mcstack
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
MpdZdcGeo.cxx
MpdZdcGeoPar.cxx
MpdZdc.cxx
MpdZdcContFact.cxx
MpdZdcDigi.cxx
MpdZdcDigiProducer.cxx
MpdZdcPoint.cxx
#MpdZdcTstSim.cxx
MpdZdcDigiPar.cxx
MpdZdcDigiScheme.cxx
)

# List of header files
Set(HEADERS
MpdZdcGeo.h
MpdZdcGeoPar.h
MpdZdc.h
MpdZdcContFact.h
MpdZdcDigi.h
MpdZdcDigiProducer.h
MpdZdcPoint.h
#MpdZdcTstSim.h
MpdZdcDigiPar.h
MpdZdcDigiScheme.h
MpdZdcPsdGeoPar.h
MpdZdcPsdDigiPar.h
)

Set(LINKDEF MpdZdcLinkDef.h)
Set(LIBRARY_NAME Zdc)
Set(DEPENDENCIES)

GENERATE_LIBRARY()
