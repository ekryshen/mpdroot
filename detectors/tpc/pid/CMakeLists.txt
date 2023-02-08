# Create a library called "libtpcPid"

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_CURRENT_LIST_DIR}
 ${CMAKE_SOURCE_DIR}/detectors/tof
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
 MpdTPCpid.cxx
 MpdParticleIdentification.cxx
)

Set(HEADERS
 MpdTPCpid.h
 MpdParticleIdentification.h # temporary workaround, this file should be somewhere else as is TPC & TOF dependent
)

Set(LINKDEF tpcPidLinkDef.h)
Set(LIBRARY_NAME tpcPid)
Set(DEPENDENCIES )

GENERATE_LIBRARY()
