set(PACKAGE TShieldGeneralGenerator) # this line comes after last add_subdirectory
# Create a library called "libTShieldGenerator"
set(INCLUDE_DIRECTORIES
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${TSHIELD_ROOT_PATH}/hadgen
  ${TSHIELD_ROOT_PATH}/shield
  ${TSHIELD_ROOT_PATH}/THadgen
  ${TSHIELD_ROOT_PATH}/TShield
)

set(SYSTEM_INCLUDE_DIRECTORIES
  ${ROOT_INCLUDE_DIR}
)

include_directories(${INCLUDE_DIRECTORIES})
Include_Directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

link_directories(
  ${ROOT_LIBRARY_DIR}
  ${FAIRROOT_LIBRARY_DIR}
)

# List of source files
set(SRCS
  TShieldGenerator.cxx
)

set(HEADERS)

set(LINKDEF ${PACKAGE}LinkDef.h)
set(LIBRARY_NAME ${PACKAGE})
set(DEPENDENCIES TShield THadgen Core Base)

GENERATE_LIBRARY()
unset(PACKAGE) # not really necessary, it will be overwritten by following CMakeLists.txt
