set(PACKAGE MpdGeneralGenerator) # this line comes after last add_subdirectory
# Create a library called "libMpdGeneralGenerator"
set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_SOURCE_DIR}/generators
  ${CMAKE_SOURCE_DIR}/gconfig
  ${CMAKE_CURRENT_SOURCE_DIR}
)

Set(SYSTEM_INCLUDE_DIRECTORIES
  ${ROOT_INCLUDE_DIR}
)

if(PLUTO_FOUND)
  set(SYSTEM_INCLUDE_DIRECTORIES ${SYSTEM_INCLUDE_DIRECTORIES} ${PLUTO_INCLUDE_DIR})
EndIf(PLUTO_FOUND)

include_directories(${INCLUDE_DIRECTORIES})
include_directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

Set(LINK_DIRECTORIES
  ${ROOT_LIBRARY_DIR}
  ${FAIRROOT_LIBRARY_DIR}
)

if(PLUTO_FOUND)
  set(LINK_DIRECTORIES ${LINK_DIRECTORIES} ${PLUTO_LIBRARY_DIR})
EndIf(PLUTO_FOUND)

link_directories(${LINK_DIRECTORIES})

set(SRCS
  MpdGeneralGenerator.cxx
)

set(HEADERS)

set(LINKDEF ${PACKAGE}LinkDef.h)
set(LIBRARY_NAME ${PACKAGE})
set(DEPENDENCIES Core Base)

GENERATE_LIBRARY()
unset(PACKAGE) # not really necessary, it will be overwritten by following CMakeLists.txt
