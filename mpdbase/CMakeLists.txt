# Create a library called "libMpdBase" which includes the source files given in the array.
# The extension is already found. Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_SOURCE_DIR}/mpdbase
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
 MpdHelix.cxx
 MpdTrack.cxx
 MpdVertex.cxx
 MpdEvent.cxx
 MpdCodeTimer.cxx
 MpdGenTrack.cxx
 MpdGenTrackTask.cxx
)

if(RULE_CHECKER_FOUND)
  CHECK_RULES("${SRCS}" "${INCLUDE_DIRECTORIES}" MPDBASE_RULES)
endif(RULE_CHECKER_FOUND)

Set(HEADERS)
Set(LINKDEF MpdBaseLinkDef.h)
Set(LIBRARY_NAME MpdBase)
Set(DEPENDENCIES )

GENERATE_LIBRARY()
