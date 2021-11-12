# The name of our project is "MPDROOT".  CMakeLists files in this project can
# refer to the root source directory of the project as ${MPDROOT_SOURCE_DIR}
# or as ${CMAKE_SOURCE_DIR} and to the root binary directory of the project as
# ${MPDROOT_BINARY_DIR} or ${CMAKE_BINARY_DIR}.

# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 3.0 FATAL_ERROR)

# Set name of our project to "MPDROOT". Has to be done
# after check of cmake version
project(MPDROOT)

IF(DEFINED ENV{SIMPATH})
  include(cmake/CMakeListsLegacy.cmake)
ELSE(DEFINED ENV{SIMPATH})
  include(cmake/CMakeLists.cmake)
ENDIF()
