# Find NLOHMANN_JSON
#
# Find the NLOHMANN_JSON include directory
# 
# if you need to add a custom library search path, do it via via CMAKE_PREFIX_PATH 
# 
# This module defines
#  NLOHMANN_JSON_INCLUDE_DIRS, where to find headers, etc.
#  NLOHMANN_JSON_FOUND, If false, do not try to use NLOHMANN_JSON.

# only look in default directories
find_path(
  NLOHMANN_JSON_INCLUDE_DIR
  NAMES "nlohmann/json.hpp"
  DOC "NLOHMANN_JSON include dir"
  HINTS ${NLOHMANN_JSON_ROOT}/include  /usr/include
)

if (NOT NLOHMANN_JSON_INCLUDE_DIR)
  MESSAGE(FATAL_ERROR "Error, NLOHMANN_JSON include dir not found. Did you set NLOHMANN_JSON_ROOT variable or have you installed json-devel package?")
endif()

set(NLOHMANN_JSON_INCLUDE_DIRS ${NLOHMANN_JSON_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set NLOHMANN_JSON_FOUND to TRUE
# if all listed variables are TRUE, hide their existence from configuration view
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NLOHMANN_JSON DEFAULT_MSG NLOHMANN_JSON_INCLUDE_DIR)
mark_as_advanced (NLOHMANN_JSON_FOUND NLOHMANN_JSON_INCLUDE_DIR)
