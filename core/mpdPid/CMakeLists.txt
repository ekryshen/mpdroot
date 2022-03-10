# Create a library called "libMpdPid" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_SOURCE_DIR}/core/mpdBase
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
  MpdPid.cxx
  MpdPidQA.cxx
)

# List of header files
set(HEADERS
  MpdPid.h
  MpdPidQA.h
  MpdPidUtils.h
)

set(LINKDEF MpdPidLinkDef.h)
set(LIBRARY_NAME MpdPid)
set(DEPENDENCIES Hist Base)

GENERATE_LIBRARY()
