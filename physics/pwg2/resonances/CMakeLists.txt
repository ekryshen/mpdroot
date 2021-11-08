# Create a library called "MpdResonances" which includes the source files given in the array.
# The extension is already found. Any number of sources could be listed here.
set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_SOURCE_DIR}/mcstack
 ${CMAKE_SOURCE_DIR}/mpdbase
 ${CMAKE_SOURCE_DIR}/mpddst
 ${CMAKE_SOURCE_DIR}/physics/pwg2/resonances
)

Set(SYSTEM_INCLUDE_DIRECTORIES
 ${ROOT_INCLUDE_DIR}
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
 MpdRoInvMassTask.cxx
)

Set(LINKDEF MpdResonancesLinkDef.h)
Set(LIBRARY_NAME MpdResonances)
Set(DEPENDENCIES Core Base MpdDst)

GENERATE_LIBRARY()