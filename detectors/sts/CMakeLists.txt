# Create a library called "libSts" which includes the source files given in
# the array.
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_SOURCE_DIR}/simulation/mcStack
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
  MpdIts5spd.cxx
  MpdItsHit5spd.cxx
  MpdItsHitProducer5spd.cxx
  MpdStsContFact.cxx
  MpdStsGeo.cxx
  MpdStsGeoPar.cxx
  MpdStsPoint.cxx
)

Set(HEADERS)
Set(LINKDEF MpdStsLinkDef.h)
Set(LIBRARY_NAME Sts)
Set(DEPENDENCIES MpdMCStack)

GENERATE_LIBRARY()
