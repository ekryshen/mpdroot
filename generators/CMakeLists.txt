# Create a library called "libMpdGen" which includes the source files given in
# the array.
# The extension is already found.  Any number of sources could be listed here.

EXEC_PROGRAM( gcc ARGS "-dumpversion" OUTPUT_VARIABLE GCC_VERSION )
STRING(REGEX REPLACE "^([0-9]+)\\.[0-9]+\\.[0-9]" "\\1" req_gcc_major_vers "${GCC_VERSION}")
STRING(REGEX REPLACE "^[0-9]+\\.([0-9]+)\\.[0-9]" "\\1" req_gcc_minor_vers "${GCC_VERSION}")
IF(${req_gcc_major_vers} MATCHES "4" AND NOT ${req_gcc_minor_vers} MATCHES "0")
ENDIF(${req_gcc_major_vers} MATCHES "4" AND NOT ${req_gcc_minor_vers} MATCHES "0")

add_subdirectory(genFactory)
add_subdirectory(unigen_format)

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_SOURCE_DIR}/config
 ${CMAKE_SOURCE_DIR}/mpdbase
 ${CMAKE_SOURCE_DIR}/mpddst/mcDst
 ${CMAKE_SOURCE_DIR}/generators
 ${CMAKE_SOURCE_DIR}/generators/unigen_format
 ${CMAKE_SOURCE_DIR}/mcstack
)

Set(SYSTEM_INCLUDE_DIRECTORIES
  ${ROOT_INCLUDE_DIR}
  ${FAIRROOT_LIBRARY_DIR}
  ${GEANT3_INCLUDE_DIR}
)

#if(PLUTO_FOUND)
#set(SYSTEM_INCLUDE_DIRECTORIES ${SYSTEM_INCLUDE_DIRECTORIES} ${PLUTO_INCLUDE_DIR})
#EndIf(PLUTO_FOUND)

include_directories(${INCLUDE_DIRECTORIES})
Include_Directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
 ${ROOT_LIBRARY_DIR}
 ${FAIRROOT_LIBRARY_DIR}
)
#if(PLUTO_FOUND)
#set(LINK_DIRECTORIES ${LINK_DIRECTORIES} ${PLUTO_LIBRARY_DIR})
#EndIf(PLUTO_FOUND)
 
link_directories(${LINK_DIRECTORIES})

# List of source files
set(SRCS
MpdMCEventHeader.cxx
MpdLAQGSMGenerator.cxx
MpdHypYPtGenerator.cxx
MpdRainGenerator.cxx
MpdPHSDGenerator.cxx
MpdUrqmdGenerator.cxx
MpdGetNumEvents.cxx
MpdHistoGenerator.cxx
Mpd3fdGenerator.cxx
MpdVHLLEGenerator.cxx
MpdFreezoutGenerator.cxx
MpdEPOSGenerator.cxx
MpdDCMSMMGenerator.cxx
MpdPHQMDGenerator.cxx
MpdMcDstGenerator.cxx
MpdDecayer.cxx
MpdTPythia8Generator.cxx
MpdSmashGenerator.cxx
MpdDecayerPyt8.cxx
MpdUnigenGenerator.cxx
)
#if(PLUTO_FOUND)
#set(SRCS ${SRCS} FairPlutoReactionGenerator.cxx MpdPlutoGenerator.cxx)
#set(SRCS ${SRCS} MpdPlutoGenerator.cxx)
#EndIf(PLUTO_FOUND)

Set(HEADERS)
Set(LINKDEF GenLinkDef.h)
Set(LIBRARY_NAME MpdGen)
Set(DEPENDENCIES Base EG MathCore Physics RIO Core MpdMcDst EGPythia6 EGPythia8 UniGenFormat)

GENERATE_LIBRARY()

