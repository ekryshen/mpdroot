# Create a library called "libUniCommon"

set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_CURRENT_LIST_DIR}
)

Set(SYSTEM_INCLUDE_DIRECTORIES
 ${Boost_INCLUDE_DIRS}
 ${ROOT_INCLUDE_DIR}
 ${XML2_INCLUDE_DIR}
)

include_directories(${INCLUDE_DIRECTORIES})
Include_Directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
 ${ROOT_LIBRARY_DIR}
 ${XML2_LIBRARY_DIR}
)

link_directories(${LINK_DIRECTORIES})

# List of source files
set(SRCS
  UniConnection.cxx
  UniValue.cxx
  UniSearchCondition.cxx
  UniGenerateClasses.cxx
  UniGeoConverter.cxx
  # add definitions of the possible database structures and enums
  db_settings.h
  db_structures.h
)

if (DEFINED LIBXML2_LIBRARIES) # if we have found XML2 library using find_package, then use this found libxml2
  SET(MY_XML2 ${LIBXML2_LIBRARIES})
else()
  SET(MY_XML2 xml2)
endif()

set(HEADERS)
Set(LINKDEF UniCommonLinkDef.h)
Set(LIBRARY_NAME UniCommon)
Set(DEPENDENCIES Core Geom Gpad Net ${MY_XML2})

GENERATE_LIBRARY()
