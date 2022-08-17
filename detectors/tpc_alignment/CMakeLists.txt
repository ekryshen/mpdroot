# Create a library called "libTpcAlignment" .

set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_SOURCE_DIR}
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

set(SRCS
  Debug.cxx
  FileHelper.cxx
  Linalg.cxx
  Alignment.cxx 
  Run.cxx 
  TpcAlignmentTask.cxx
)

Set(HEADERS)

Set(LINKDEF tpcAlignmentLinkDef.h)
Set(LIBRARY_NAME TpcAlignment)
Set(DEPENDENCIES Core Base )

GENERATE_LIBRARY()