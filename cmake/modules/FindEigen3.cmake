# Find Eigen3
#
# Find the Eigen3 include directory
# 
# if you need to add a custom library search path, do it via via CMAKE_PREFIX_PATH 
# 
# This module defines
#  Eigen3_INCLUDE_DIRS, where to find headers, etc.
#  Eigen3_FOUND, If false, do not try to use Eigen3.

# only look in default directories
find_path(
	Eigen3_INCLUDE_DIR
	NAMES "signature_of_eigen3_matrix_library"
	DOC "FMT include dir"
	HINTS ${EIGEN3_ROOT}/include/eigen3  /usr/include/eigen3
)

if (NOT Eigen3_INCLUDE_DIR)
  MESSAGE(FATAL_ERROR "Error, Eigen3 include dir not found. Did you set EIGEN3_ROOT variable?")
endif()

set(Eigen3_INCLUDE_DIRS ${Eigen3_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set Eigen3_FOUND to TRUE
# if all listed variables are TRUE, hide their existence from configuration view
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eigen3 DEFAULT_MSG Eigen3_INCLUDE_DIR)
mark_as_advanced (Eigen3_FOUND Eigen3_INCLUDE_DIR)
