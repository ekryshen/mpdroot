# Create a library called "libmulti"

set(INCLUDE_DIRECTORIES
 ${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_CURRENT_LIST_DIR}
 ${CMAKE_SOURCE_DIR}/detectors/tof
 ${CMAKE_SOURCE_DIR}/detectors/tpc/pid
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
 MpdParticleIdentification.cxx
)

Set(HEADERS
 MpdParticleIdentification.h 
)

Set(LINKDEF multiLinkDef.h)
Set(LIBRARY_NAME multi)
Set(DEPENDENCIES tpcPid)

GENERATE_LIBRARY()
