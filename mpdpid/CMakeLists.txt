# Create a library called "libMpdPid" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_SOURCE_DIR}/mpdbase
 ${CMAKE_SOURCE_DIR}/mpdpid
)

Set(SYSTEM_INCLUDE_DIRECTORIES
 ${ROOT_INCLUDE_DIR}
)

include_directories(${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
 ${ROOT_LIBRARY_DIR}
 ${FAIRROOT_LIBRARY_DIR}
)
 
link_directories(${LINK_DIRECTORIES} )

# List of source files
set(SRCS
MpdPid.cxx
MpdPidQA.cxx
)

# List of header files
Set(HEADERS
MpdPidUtils.h
MpdPid.h
MpdPidQA.h
)
Set(LINKDEF MpdPidLinkDef.h)
Set(LIBRARY_NAME MpdPid)
Set(DEPENDENCIES)

GENERATE_LIBRARY()
