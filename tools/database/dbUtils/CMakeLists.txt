# Create a library called "libDbUtils"

set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_SOURCE_DIR}/tools/database/uniCommon
  ${CMAKE_SOURCE_DIR}/tools/database/uniDb
)

Set(SYSTEM_INCLUDE_DIRECTORIES
  ${Boost_INCLUDE_DIRS}
  ${ROOT_INCLUDE_DIR}
  ${XML2_INCLUDE_DIR}
)

include_directories(${INCLUDE_DIRECTORIES})
include_directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
  ${ROOT_LIBRARY_DIR}
  ${XML2_LIBRARY_DIR}
)

link_directories(${LINK_DIRECTORIES})

if (DEFINED LIBXML2_LIBRARIES) # if we have found XML2 library using find_package, then use this found libxml2
  SET(MY_XML2 ${LIBXML2_LIBRARIES})
else()
  SET(MY_XML2 xml2)
endif()

# List of source files
set(SRCS
  DbUtils.cxx
)

Set(HEADERS
  DbUtils.h
)
Set(LINKDEF dbUtilsLinkDef.h)
Set(LIBRARY_NAME DbUtils)
Set(DEPENDENCIES UniCommon UniDb Core Geom Gpad Net ${MY_XML2})

GENERATE_LIBRARY()

# add_executable (show_simulation_files show_simulation_files.cxx)
# target_link_libraries(show_simulation_files ${ROOT_LIBRARIES})

# add_executable (show_experiment_files show_experiment_files.cxx)
# target_link_libraries(show_experiment_files ${ROOT_LIBRARIES})
