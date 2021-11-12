# Find CLHEP
#
# Find the CLHEP includes
# 
# if you nee to add a custom library search path, do it via via CMAKE_PREFIX_PATH 
# 
# This module defines
#  CLHEP_INCLUDE_DIRS, where to find header, etc.
#  CLHEP_FOUND, If false, do not try to use CLHEP.

# only look in default directories
find_path(
	CLHEP_INCLUDE_DIR
	NAMES CLHEP/ClhepVersion.h
	DOC "CLHEP include dir"
	HINTS ${CLHEP_ROOT}/include
)

if (NOT CLHEP_INCLUDE_DIR)
  MESSAGE(FATAL_ERROR "Error, CLHEP include dir not found. Did you set CLHEP_ROOT variable?")
endif()

set(CLHEP_INCLUDE_DIRS ${CLHEP_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set CLHEP_FOUND to TRUE
# if all listed variables are TRUE, hide their existence from configuration view
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CLHEP DEFAULT_MSG CLHEP_INCLUDE_DIR)
mark_as_advanced (CLHEP_FOUND CLHEP_INCLUDE_DIR)
