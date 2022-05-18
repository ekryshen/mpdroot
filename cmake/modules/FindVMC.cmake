# Find VMC
#
# Find the VMC includes
# 
# if you nee to add a custom library search path, do it via via CMAKE_PREFIX_PATH 
# 
# This module defines
#  VMC_INCLUDE_DIRS, where to find header, etc.
#  VMC_LIB_DIR, where to find library
#  VMC_FOUND, If false, do not try to use VMC.

# only look in default directories
find_path(
  VMC_INCLUDE_DIR
  NAMES vmc/TMCtls.h
  DOC "VMC include dir"
  HINTS ${VMC_ROOT}/include
)

find_path(VMC_LIB_DIR NAMES libVMC.so PATHS
  ${VMC_ROOT}/lib
  NO_DEFAULT_PATH
)

if (NOT VMC_INCLUDE_DIR)
  MESSAGE(WARNING "Warning, VMC include dir not found. Did you set VMC_ROOT variable?")
endif()

set(VMC_INCLUDE_DIRS ${VMC_INCLUDE_DIR} ${VMC_INCLUDE_DIR}/vmc)

# handle the QUIETLY and REQUIRED arguments and set VMC_FOUND to TRUE
# if all listed variables are TRUE, hide their existence from configuration view
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(VMC DEFAULT_MSG VMC_INCLUDE_DIR VMC_LIB_DIR)
mark_as_advanced (VMC_FOUND VMC_INCLUDE_DIR VMC_LIB_DIR)
