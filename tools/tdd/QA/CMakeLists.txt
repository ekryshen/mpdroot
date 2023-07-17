# Create a library called "libQA" which is library for generating QA plots

set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_SOURCE_DIR}/detectors/tpc
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
 BaseQA.cpp
 QA_TpcClusterHitFinder.cxx
)

Set(HEADERS
 BaseQA.h
 QA_TpcClusterHitFinder.h
)


Set(LINKDEF qaLinkDef.h)
Set(LIBRARY_NAME QA)
Set(DEPENDENCIES )

GENERATE_LIBRARY()
