set(INCLUDE_DIRECTORIES
${BASE_INCLUDE_DIRECTORIES}
${CMAKE_SOURCE_DIR}/mpddst/mcDst
)

Set(SYSTEM_INCLUDE_DIRECTORIES
 ${ROOT_INCLUDE_DIR}
 ${GSL_INCLUDE_DIRS}
 ${CLHEP_INCLUDE_DIR}
 ${Boost_INCLUDE_DIRS}
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
MpdMcRun.cxx
MpdMcPIDConverter.cxx
MpdMcParticle.cxx
MpdMcEvent.cxx
MpdMcDstReader.cxx
MpdMcDst.cxx
MpdMcArrays.cxx
)

Set(HEADERS)
Set(LINKDEF MpdMcDstLinkDef.h)
Set(LIBRARY_NAME MpdMcDst)
Set(DEPENDENCIES Core Base)

GENERATE_LIBRARY()
