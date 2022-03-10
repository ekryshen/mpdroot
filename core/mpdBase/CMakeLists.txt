# Create a library called "libMpdBase" which includes the source files given in the array.
# The extension is already found. Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_CURRENT_LIST_DIR}
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
  MpdCodeTimer.cxx
  MpdEvent.cxx
  MpdGenTrack.cxx
  MpdGenTrackTask.cxx
  MpdHelix.cxx
  MpdTrack.cxx
  MpdVertex.cxx
)

if(RULE_CHECKER_FOUND)
  CHECK_RULES("${SRCS}" "${INCLUDE_DIRECTORIES}" MPDBASE_RULES)
endif(RULE_CHECKER_FOUND)

set(HEADERS
  MpdBase.h
  MpdCodeTimer.h
  MpdEvent.h
  MpdGenTrack.h
  MpdGenTrackTask.h
  MpdHelix.h
  MpdTrack.h
  MpdVertex.h
)

set(LINKDEF MpdBaseLinkDef.h)
set(LIBRARY_NAME MpdBase)
set(DEPENDENCIES)

GENERATE_LIBRARY()
