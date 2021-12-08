set(INCLUDE_DIRECTORIES
${BASE_INCLUDE_DIRECTORIES}
 ${CMAKE_SOURCE_DIR}/generators
 ${CMAKE_SOURCE_DIR}/generators/genFactory
 ${CMAKE_SOURCE_DIR}/shield_pack/THadgen
 ${CMAKE_SOURCE_DIR}/shield_pack/hadgen
 ${CMAKE_SOURCE_DIR}/macro/mpd # needed to load loadlibs.C
 ${CMAKE_SOURCE_DIR} # needed to find gconfig/basiclibs.C
)

Set(SYSTEM_INCLUDE_DIRECTORIES
 ${ROOT_INCLUDE_DIR}
 ${FAIRROOT_LIBRARY_DIR}
)

include_directories(${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
 ${ROOT_LIBRARY_DIR}
 ${FAIRROOT_LIBRARY_DIR}
)
 
link_directories(${LINK_DIRECTORIES})

# List of source files
set(SRCS
MpdGeneratorsFactory.cxx
)

Set(HEADERS)
Set(LINKDEF LinkDefFactory.h)
Set(LIBRARY_NAME MpdGenFactory)
Set(DEPENDENCIES)

GENERATE_LIBRARY()
