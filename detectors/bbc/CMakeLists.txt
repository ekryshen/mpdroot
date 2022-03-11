# Create a library called "libBbc" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
  ${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_SOURCE_DIR}/simulation/mcStack
)

Set(SYSTEM_INCLUDE_DIRECTORIES
  ${ROOT_INCLUDE_DIR}
)

include_directories(${INCLUDE_DIRECTORIES})
include_Directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
  ${ROOT_LIBRARY_DIR}
)

link_directories(${LINK_DIRECTORIES})

# List of source files
set(SRCS
  MpdBbc.cxx
  MpdBbcContFact.cxx
  MpdBbcGeo.cxx
  MpdBbcGeoPar.cxx
  MpdBbcHit.cxx
  MpdBbcHitProducer.cxx
  MpdBbcPoint.cxx
)

Set(HEADERS)
Set(LINKDEF MpdBbcLinkDef.h)
Set(LIBRARY_NAME Bbc)
Set(DEPENDENCIES Core Base MpdMCStack)

GENERATE_LIBRARY()
